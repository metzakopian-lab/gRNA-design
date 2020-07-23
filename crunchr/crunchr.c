#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>

enum {
  NBIT_PER_BYTE = 8,
  NBIT_HASH = 20, /* number of hashed bits */
  MEMBLKSZ = 1200,
  LINEWIDTH_MAX = 1024, /* maximum line width */
};

enum ERROR_CODES {
    ERRCODE_SUCCESS = 0,
    ERRCODE_FAILURE = -1,
    ERRCODE_NOMEM = 2,     /* memory allocation failed */
    ERRCODE_SEQLEN = 3,    /* wrong sequence legth */
};


/** yoonha:
 * stagger TGACTGA
 * TGACTGATCTTGTGGAAAGGACGAAACACCGNNNNNNNNNNNNNNNNNNNN
 * TGACTGATCTTGTGCAAAGGACGAAACACCGCAAGCCGCTCATCGTGCAT
 * TGACTGATCTTGTGGAAAGGACGCAACACCGTCTTTCTCGTCGTCGAATG
 * stagger ACTGA
 * ACTGATCTTGTGGAAAGGACGAAACACCGNNNNNNNNNNNNNNNNNNNNG 
 * ACTGATCTTGTGGAATAGACGAAACAGCGTACACAACTGCTTCCAGGTAG
 * ACTGATCGTGTGGAAGGGACGAAAGACCGGATTATCGGGTGTGTGCTGAG
 *
 * stagger GA
 * GATCTTGTGGAAAGGACGAAACACCGNNNNNNNNNNNNNNNNNNNNGTTT
 * GATCTTGTGGAAAGGACGAAACACCGAGGTTCTCTTGACCCCCTGGGTTT
 * GATCTTGTGGAAAGCACGAAACACCGTGCCAAGTTCCAACTTGTTAGTAT
 *
 * stagger None
 * TCTTGTGGAAAGGACGAAACACCGNNNNNNNNNNNNNNNNNNNNGTTTTA
 * TCTTGTGGAAAGGACGAAACACCGACGCCGTCGATCATATCTGCGTTTTA
 * TCTTGTGGAAAGGACGAAACACCGTCTCTAGAACTGCACCATTCGTTTTA
 * TCTTGTGGAAAGGACGAAACACCGCAGCATGTTTCACGTTCGAAGTTTTA
 *
 * 
 *
 * GATCTTGTGGAAAGGACGAAACACCGAGGTTCTCTTGACCCCCTGGGTTT
 *
 * GCTCGTCGAAAGGACGAAACACCGGCGCAAAACTATTCCCATACGTTTTA
 * GCTCGTCGAAAGGACACGACACCCCCCTGCGTATCAGCTTCGACGCTATA
 * TCTTATGTCGAGGAAAAAACACCCAGATATGGAAGAGCTTGTTGTTGTCG
 * TCTTTTGGCAAGGACGAAACACCACCCTCGCACGCAGACGGCCAGTCTTC
 * TATTTTGTCAAGTACTAAACACCAGCCACTTGCTCTGCTTCGCCTTCTTC
 * GATCTTGTGGAAAGGACGAAACACCGAGGTTCTCTTGACCCCCTGGGTTT
 */

typedef struct HASHELEM {
  uint32_t nhk;
  uint64_t count;
  uint32_t idx_next;
} HASHELEM;

typedef struct GuideHash {
  HASHELEM *elemp;
  uint32_t n_elem;
  uint32_t n_alloc;
  uint32_t ix_free;
} GuideHash;

const char *errMsgString(int errcode) 
{
  switch(errcode) {
  case ERRCODE_SUCCESS:
    return "success";
  case ERRCODE_FAILURE:
    return "failure";
  case ERRCODE_NOMEM:
    return "memory allocation failed";
  case ERRCODE_SEQLEN:
    return "wrong sequence length";
  default:
    return "unknown error code";
  }
}

int errMsgPrint(int errc, const char *fnam, int lineno)
{
  return fprintf(stderr, "ERROR: %s in file %s, line %i", errMsgString(errc), fnam, lineno);
}

#define ERRMSG_PRINT(errc) errMsgPrint((msg), __FILE__, __LINE__)

void destroyGuideHash(GuideHash *ghp)
{
  if (NULL != ghp) {
    free(ghp->elemp);
    free(ghp);
  }
}

GuideHash *setupGuideHash(void)
/**< Setup a hash index for a sequence length of ntnum nucleotides 
 */
{
  GuideHash *ghp = NULL;
 
  if (NULL != (ghp = malloc(sizeof(GuideHash)))) {    
    ghp->n_elem = (((uint32_t) 1) << NBIT_HASH) - 1;
    ghp->ix_free = ghp->n_elem;
    ghp->n_alloc = ghp->n_elem<<1;
    ghp->elemp = calloc(ghp->n_alloc, sizeof(HASHELEM));
    if (NULL == ghp->elemp) {
      destroyGuideHash(ghp);
      ghp = NULL;
    }
  }
  return ghp;
}

int setGuideHashKey(GuideHash * const ghp, const char *keyp)
{
  ;
}

GuideHash *loadGuideHashFromFile(const char *fnam)
/**< load the guide sequences from a *.csv file, tab delimited,
 * with the guide sequence in the first column 
 */
{
  char linbuf[LINEWIDTH_MAX];
  GuideHash *ghp = NULL;
  FILE *fp = fopen(fnam, "r");

  if (NULL == fp) {
    fprintf(stderr, "Could not open file '%s'", fnam);
    return ghp;
  }

  ghp = setupGuideHash();
  if (NULL == ghp) 
    return ghp;
  
  is_full_line = 1;
  while (NULL != fgets(linbuf, LINEWIDTH_MAX, fp)) {
    size_t nc = strlen(linbuf);
    if (is_full_line) {
      char *fldp = strtok(linbuf, '\t');
      if (*fldp != '#') {
	setGuideHashKey(ghp, fldp);
      }
    }
    is_full_line = nc > 0 && linbuf[nc-1] == '\n';
  }
  if (ferror(fp) != 0) {
    perror("error when reading file.");
    destroyGuideHash(ghp);
    ghp = NULL;
  }

  return ghp;
}

int main(int argc, char *argv[]) {
  char *infnamp = NULL;
  samFile *sfp; /* samFile is an alias for htsFile */
  bam_hdr_t *headp;
  bam1_t *bamp;
  int8_t *sbufp = NULL;
  uint32_t sbuflen = 0;

  if (argc != 2) {
    fprintf(stderr, "usage: %s <input file [BAM|SAM|CRAM]>\n", argv[0]);
    return 1;
  }

  infnamp= argv[1];
  sfp = sam_open(infnamp, "r"); /* alias for hts_open */
  if (NULL == sfp) {
    fprintf(stderr, "Cannot read file \"%s\"", infnamp);
    return 1;
  }

  headp = sam_hdr_read(sfp); /* declared in htslib/sam.h */
  bamp = bam_init1();        /* declared in htslib/sam.h */

  while (sam_read1(sfp, headp, bamp) >= 0) { /* sam_read1 declared in htslib/sam.h */
    if (!(bamp->core.flag & BAM_FREAD2)) {
      int32_t i, qlen = bamp->core.l_qseq;
      uint8_t* seqp = bam_get_seq(bamp);
      if (sbuflen <= (uint32_t) qlen) {
	if (NULL != sbufp)
	  free(sbufp);
	sbuflen = (qlen + MEMBLKSZ + 1)/MEMBLKSZ;
	sbuflen *= MEMBLKSZ;
	sbufp = malloc(sbuflen);
      }
      for (i = 0; i < qlen; ++i)
	sbufp[i] = seq_nt16_str[bam_seqi(seqp, i)];
      sbufp[i] = '\0';
      printf("%s %i\n", (char *) sbufp, qlen);
    }
  }

  bam_destroy1(bamp);
  bam_hdr_destroy(headp);

  sam_close(sfp); /* alias for hts_close */
  
  return 0;
}
