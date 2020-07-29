/** Design a guide library */
//#define DESIGUIDE_DEBUG
#define DESIGUIDE_VERBOSE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>

#include "smalt/elib.h"
#include "smalt/array.h"
#include "smalt/sequence.h"


enum {
  PAMSLEN = 3,
  SEEDCORE = 18, /* number of bp upsteam of pam minimum number of mismatches */
  SIZE_STANDARD_ALPHABET = 4,
  NBITS_KEYCOD = 2,
  NBITS_ALPHABET = 3,
  MAXN_PER_UNIT = 10,
  NBASES_WORD = 20,
  MEMALLOC_BLKSZ = 8192,
  FILHEAD_NUMFLD = 2,
  CHRNAM_MAX = 64,
  GENID_MAX = 64,
  LINBUF_MAX = 512,
  STRAND_FORWARD = 1,
  STRAND_REVERSE_COMPLEMENT = -1, 
  OFFTARGETS_MAXNUM_REPORT = 10, /* maximum number of off-targets to report */
  MINIMUM_EDIT_DISTANCE = 1,
  NAMBUFSZ = 64,
  FILELIN_MAX = 256
};

static const char const DESIGUIDE_VERSION[] = 
  "$Id: desiguide.c,v 1.5 2017/03/28 11:03:44 hp3 Exp $";

typedef uint32_t GUIDEPOS_T;
typedef uint64_t GUIDEWORD_T;
typedef uint64_t HEADER_T;
typedef unsigned long long int ULL_T;

typedef struct GuideIndex_ {
  GUIDEWORD_T *guidr;  /**< Array of 2-bit encoded guide sequences (excl PAM) */
  GUIDEPOS_T *posr;    /**< Array of positions, 0-based in sequence set
			* position is 1st base of word along reference strand */
  char pam[PAMSLEN+1];
} GuideIndex;

typedef struct BED_ {
  char chrnam[CHRNAM_MAX]; /**< chromosome name */
  SEQLEN_t start;          /**< 1-based start position within chromosome */
  SEQLEN_t end;            /**< 1-based end position within chromosome */
  int strand;              /**< 1: forward, -1: reverse complement */
  SETSIZ_t os;             /**< 0-based offset of 1st base of segment in set */
  SETSIZ_t oe;             /**< 0-based offset of last base of segment in set */
  SEQNUM_t seqidx;         /**< sequence index in set */
  char genID[GENID_MAX];  /**< gene ID */
  char uniprotID[GENID_MAX];
  char ensemblID[GENID_MAX];
  char annotation[FILELIN_MAX];
} BED;

typedef struct MATCH_ {
  uint32_t idx;    /* index of matching word */
  GUIDEWORD_T df;  /* bitstring 2 bits per base == 00 match, != 00 mismatch */
  GUIDEWORD_T gw;  /* actual guide */
  int nmm;         /* number of mismatches */
} MATCH;

static const GUIDEWORD_T GUIDE_MASK = 
  (((GUIDEWORD_T) 1)<<NBASES_WORD*NBITS_KEYCOD) - 1;
static const GUIDEWORD_T GUIDE_RCFLG = ((GUIDEWORD_T) 1)<<63;
static const uint32_t ALPHABET_MASK = 0x07;
static const char CODEC_ALPHABET[] = "ACGTXN";
static const uint32_t FILHEAD_SIGNATURE = 0xFE7310;
static const size_t GUIDEPOS_MAX = UINT32_MAX;
static const char FILEXT_GUIDES[] = ".gid";
static const char FILEXT_POS[] = ".pos";
static const char NOSTR[] = "NONE";

#define REVERSE_COMPLEMENT_WORD(word, nblen)				\
  do {									\
    GUIDEWORD_T w = ~word;						\
    int i = nblen;							\
    word = 0LL;								\
    while (i-- > 0) {							\
      word = (word<<NBITS_KEYCOD) | (w & SEQCOD_STDNT_MASK);		\
      w >>= NBITS_KEYCOD;						\
    }									\
  } while (0);

static char * getShortSequenceName(char *buffer, int bufsz, const char *name)
{
  char *p = buffer;
  for (; bufsz>1 && !isspace(*name) && '\0' != *name; bufsz--) {
    *p++ = *name++;
  }
  *p = '\0';

  return buffer;
}

static char * decodeWord(char bufp[], GUIDEWORD_T word, int len)
{
  int nsu;
  char *cp;
  for (nsu = len-1, cp = bufp; nsu >= 0; nsu--, cp++)
    *cp = CODEC_ALPHABET[(word>>nsu*NBITS_KEYCOD)&SEQCOD_STDNT_MASK];
  *cp = '\0';
  return bufp;
}

static int encodeWord(GUIDEWORD_T * const word,
		      const SeqFastq * const sqp)
{
  char cod;
  SEQLEN_t slen;
  int sl, nsu;
  
  const uint32_t *basp = (uint32_t *) seqFastqGetConstSequence(sqp, &slen, &cod);

  if (cod != SEQCOD_COMPRESSED)
    return ERRCODE_SEQCODE;

  sl = (slen > NBASES_WORD)? NBASES_WORD: (int) slen;
  
  *word = 0LL;
  for (nsu = MAXN_PER_UNIT-1; sl > 0; sl--) {
    uint32_t bc = ((*basp)>>(nsu*NBITS_ALPHABET))&SEQCOD_ALPHA_MASK;
    if (--nsu < 0) {
      basp++;
      nsu = MAXN_PER_UNIT-1;
    }
    *word = ((*word)<<NBITS_KEYCOD)|(bc&SEQCOD_STDNT_MASK);
  }

  *word &= GUIDE_MASK;

  return ERRCODE_SUCCESS;
}

static int scoreMatch(char * const mmstr, const MATCH * const mp)
{
  int p, sc = 0;
  GUIDEWORD_T df;

  for (p = NBASES_WORD, df = mp->df; df != 0LL && p>0;p--) {
    if ((df&SEQCOD_STDNT_MASK)) {
      sc += p;
      if (NULL != mmstr)
	mmstr[p-1] = 'X';
    } else if (NULL != mmstr)
      mmstr[p-1] = '-';
    df >>= NBITS_KEYCOD;
  }
  if (NULL != mmstr) {
    while (p>0)
      mmstr[--p] = '-';
    mmstr[NBASES_WORD] = '\0';
  }
  return sc;
}

static char *convertSeq2RY(char *sqp)
{
  char *s;
  for (s = sqp;*s != '\0';s++) {
    if ('A' == *s || 'a' == *s ||
	'G'  == *s || 'g' == *s) {
      *s = 'R';
    } else if ('C' == *s || 'c' == *s ||
	       'T'  == *s || 't' == *s) {
      *s = 'Y';
    }
  }
  return sqp;
}

static int generateMisMatchString(char * const sqp,
				  int * const mmscore,
				  const char * const s1p,
				  const char * const s2p)
{
  int i;
  int nmm = 0;
  int mms = 0;
  for (i=0; '\0' != s1p[i] && '\0' != s2p[i]; i++) {
    if (s1p[i] == s2p[i]) {
      sqp[i] = '-';
    } else if ('N' == s1p[i] || 'N' == s2p[i]) {
      sqp[i] = 'N';
    } else {
      sqp[i] = 'X';
      nmm++;
      mms += i+1;
    }
  }
  sqp[i] = '\0';
  if (mmscore) {
    *mmscore = mms;
  }

  return nmm;
}

static int cmpBED(const void *ap, const void *bp) {
  int rv;
  const BED *a = (const BED *) ap;
  const BED *b = (const BED *) bp;
  if (a->os < b->os)
    rv = -1;
  else if (a->os > b->os)
    rv = 1;
  else if (a->oe > b->oe)
    rv = -1;
  else if (a->oe < b->oe )
    rv = 1;
  else
    rv = 0;

  return rv;
}

static int initGuideIndex(GuideIndex * const gxp, 
			  size_t blksz)
{
  int errcode = ERRCODE_SUCCESS;

  if (NULL == ARRCREATE(gxp->guidr, blksz) ||
      NULL == ARRCREATE(gxp->posr, blksz)) {
    errcode = ERRCODE_NOMEM;
  } 
  return errcode;
}

static void resetGuideIndex(const GuideIndex * const gxp)
{
  ARRLEN(gxp->guidr) = 0LL;
  ARRLEN(gxp->posr) = 0LL;
}
static char *setGuideIndexPAM(GuideIndex * const gxp, 
			      const char * const pamseqp)
{
  strncpy(gxp->pam, pamseqp, PAMSLEN);
  gxp->pam[PAMSLEN] = '\0';
  return gxp->pam;
}

static void freeGuideIndex(GuideIndex * const gxp)
{
  if (gxp != NULL) {
    ARRDELETE(gxp->guidr);
    ARRDELETE(gxp->posr);
    gxp->guidr = NULL;
    gxp->posr = NULL;
  }
}

static int writeGuideIndexToFile(const GuideIndex * const gxp,
				 int writePositions,
				 const char * const fnamp)
{
  char fnambuf[FILENAME_MAX];
  HEADER_T header[FILHEAD_NUMFLD];
  HEADER_T sigfld;
  char *pamseqp = (char *) (((uint32_t *) header)+1);
  FILE *fp = NULL;
  const char *fnamextp = (writePositions)? FILEXT_POS: FILEXT_GUIDES;
  size_t fnamextlen = strlen(fnamextp);


  strncpy(fnambuf, fnamp, FILENAME_MAX - fnamextlen - 1);
  fnambuf[FILENAME_MAX - fnamextlen - 1] = '\0';
  strcat(fnambuf, fnamextp);

  if (!(fp = EFOPEN(fnambuf, "wb")))
    return ERRCODE_NOFILE;
  
  strncpy(pamseqp, gxp->pam, PAMSLEN);
  pamseqp[PAMSLEN] = '\0';
  sigfld = FILHEAD_SIGNATURE;
  if ((writePositions)) 
    sigfld &= 0x01;
  header[0] = sigfld<<32;
  strncpy(pamseqp, gxp->pam, PAMSLEN);
  pamseqp[PAMSLEN] = '\0';
  header[1] = ARRLEN(gxp->guidr);
  fwrite(header, sizeof(header[0]), FILHEAD_NUMFLD, fp);
  if ((writePositions)) {
    fwrite(gxp->posr, sizeof(GUIDEPOS_T), header[1], fp);
  } else {
    fwrite(gxp->guidr, sizeof(GUIDEWORD_T), header[1], fp);
  }

  return EFCLOSE(fp);
}

static int loadGuideIndexFromFile(GuideIndex * const gxp,
				  int loadPositions,
				  const char * const fnamp)
{
  char fnambuf[FILENAME_MAX];
  HEADER_T header[FILHEAD_NUMFLD];
  HEADER_T sigfld;
  char *pamseqp = (char *) (((uint32_t *) header)+1);
  size_t n_guides;
  void *hp;
  FILE *fp = NULL;
  const char *fnamextp = (loadPositions)? FILEXT_POS: FILEXT_GUIDES;
  size_t fnamextlen = strlen(fnamextp);

  strncpy(fnambuf, fnamp, FILENAME_MAX - fnamextlen - 1);
  fnambuf[FILENAME_MAX - fnamextlen - 1] = '\0';
  strcat(fnambuf, fnamextp);

  if (!(fp = EFOPEN(fnambuf, "rb")))
    return ERRCODE_NOFILE;

  if (fread(header, sizeof(header[0]), FILHEAD_NUMFLD, fp) != FILHEAD_NUMFLD)
    return ERRCODE_READERR;

  sigfld = header[0]>>32;
  if ((!(loadPositions) && sigfld != FILHEAD_SIGNATURE) &&
      ((loadPositions) && sigfld != FILHEAD_SIGNATURE+1))
    return ERRCODE_FILTYP;

  setGuideIndexPAM(gxp, pamseqp);
  
  n_guides = header[1];

  if (loadPositions) {
    if (NULL == (hp = ARREALLOC(gxp->posr, n_guides)))
      return ERRCODE_NOMEM;
    gxp->posr = hp;
    if (fread(gxp->posr, sizeof(GUIDEPOS_T), n_guides, fp) != n_guides)
      return ERRCODE_READERR;
    ARRLEN(gxp->posr) = n_guides;
  } else {
    if (NULL == (hp = ARREALLOC(gxp->guidr, n_guides)))
      return ERRCODE_NOMEM;
    gxp->guidr = hp;
    if (fread(gxp->guidr, sizeof(GUIDEWORD_T), n_guides, fp) != n_guides)
      return ERRCODE_READERR;
    ARRLEN(gxp->guidr) = n_guides;
  }

  return EFCLOSE(fp);
}

#ifdef DESIGUIDE_DEBUG
static int checkGuideIndex(const GuideIndex * const gxp,
			   const SeqSet * const ssp)
{
  int errcode = ERRCODE_SUCCESS;
  size_t g = ARRLEN(gxp->guidr);
  char cbuf[NBASES_WORD + 1];
  SeqFastq *sqp = seqFastqCreate(0, SEQTYP_FASTA);
  SeqCodec *codecp = seqCodecCreate();

  if (ARRLEN(gxp->posr) != g)
    return ERRCODE_FAILURE;

  while (g-- > 0) {
    int i;
    const char *seqstrp;
    SEQLEN_t slen;
    GUIDEWORD_T word = gxp->guidr[g];
    GUIDEPOS_T pos = gxp->posr[g];
    SETSIZ_t start = (word & GUIDE_RCFLG)? pos - PAMSLEN: pos;
    SETSIZ_t end = start + NBASES_WORD + PAMSLEN - 1;

    if ((errcode = seqSetFetchSegment(sqp, &start, &end, ssp, codecp)))
      return errcode;

    if ((word & GUIDE_RCFLG)) {
      if ((errcode = seqFastqReverse(sqp, codecp)))
	return errcode;
    }
    seqstrp = seqFastqGetConstSequence(sqp, &slen, NULL);
    if (slen != NBASES_WORD + PAMSLEN)
      return ERRCODE_FAILURE;
    decodeWord(cbuf, word, NBASES_WORD);
    /* check guide */
    for (i=0; i<NBASES_WORD && cbuf[i] == seqstrp[i]; i++);
    if (i < NBASES_WORD)
      return ERRCODE_FAILURE;
    seqstrp += i;
    for (i=0; 
	 i<PAMSLEN && (seqstrp[i] == gxp->pam[i] || gxp->pam[i] == 'N'); 
	 i++);
    if (i < PAMSLEN)
      return ERRCODE_FAILURE;
  }

  seqFastqDelete(sqp);
  seqCodecDelete(codecp);

  return errcode;
}

static int checkGuideIndexComplete(const GuideIndex * const gxp,
				   const SeqSet * const ssp)
{
  int errcode = ERRCODE_SUCCESS;
  int nsu;
  char cod;
  unsigned char ssflg;
  SETSIZ_t ctr, siz;
  uint32_t const *basep = (const uint32_t *) seqSetGetBaseData(&ssflg, &cod, &siz, ssp);
  size_t s = 0;
  const size_t n_guides = ARRLEN(gxp->guidr);
  GUIDEWORD_T guide = gxp->guidr[0], word = 0LL;
  SeqFastq *sqp = seqFastqCreate(0, SEQTYP_FASTA);
  SeqCodec *codecp = seqCodecCreate();

  if (cod != SEQCOD_COMPRESSED || !(ssflg&SEQSET_COMPRESSED))
    return ERRCODE_SEQCODE;
  
  if (guide & GUIDE_RCFLG) {
    REVERSE_COMPLEMENT_WORD(guide, NBASES_WORD);
    guide |= GUIDE_RCFLG;
  }

  for (nsu = MAXN_PER_UNIT-1, ctr = 0LL; ctr<siz; ctr++) {
    uint32_t bc = ((*basep)>>(nsu*NBITS_ALPHABET))&SEQCOD_ALPHA_MASK;
    if (--nsu < 0) {
      basep++;
      nsu = MAXN_PER_UNIT-1;
    }
    word = (word<<NBITS_KEYCOD)|(bc&SEQCOD_STDNT_MASK);
    if (0LL == ((word^guide)&GUIDE_MASK)) {
      /* found guide, check guide, check PAM */
      int i;
      char cbuf[NBASES_WORD + 1];
      const char *seqstrp;
      SEQLEN_t slen;
      SETSIZ_t os = ctr - NBASES_WORD + 1, oe = ctr;
      /* check guide */
      if ((errcode = seqSetFetchSegment(sqp, &os, &oe, ssp, codecp)))
	return errcode;
      seqstrp = seqFastqGetConstSequence(sqp, &slen, NULL);
      decodeWord(cbuf, word, NBASES_WORD);
      for (i=0; i<NBASES_WORD && cbuf[i] == seqstrp[i]; i++);
      if (i < NBASES_WORD && slen == NBASES_WORD) {
	/* not a match, check whether there are non-standard NTs */
	for (i=0; i<NBASES_WORD ; i++) {
	  int j;
	  for(j=0; 
	      j<SIZE_STANDARD_ALPHABET && CODEC_ALPHABET[j] != seqstrp[i]; 
	      j++);
	  if (j >= SIZE_STANDARD_ALPHABET)
	    break;
	}
      }
      if (i < NBASES_WORD)
	continue;
      /* check PAM */
      if (guide&GUIDE_RCFLG) {
	os = ctr - NBASES_WORD - PAMSLEN + 1;
      } else {
	os = ctr + 1;
      }
      oe = os + PAMSLEN - 1;
      
      if ((errcode = seqSetFetchSegment(sqp, &os, &oe, ssp, codecp)))
	return errcode;

      if (os + PAMSLEN - 1 != oe)
	continue;
      
      if ((guide&GUIDE_RCFLG) &&
	  (errcode = seqFastqReverse(sqp, codecp)))
	return errcode;

      seqstrp = seqFastqGetConstSequence(sqp, &slen, NULL);
      if (slen != PAMSLEN) 
	continue;
      
      for (i=0; i<PAMSLEN && seqstrp[i] == gxp->pam[i]; i++);
      if (i >= PAMSLEN) {
	/* PAM matches */
	if (++s < n_guides) {
	  guide = gxp->guidr[s];
	  if (guide & GUIDE_RCFLG) {
	    REVERSE_COMPLEMENT_WORD(guide, NBASES_WORD);
	    guide |= GUIDE_RCFLG;
	  }
	} else if (s == n_guides) {
	  guide = 0LL;
	} else {
	  return ERRCODE_FAILURE;
	}	
      }
    }
#ifdef DESIGUIDE_VERBOSE
    if ( 0 == ctr % 0x07FFFF)
      printf("%5.2f %% processed ...\n", ((double) ctr)*100/siz);
#endif
  }
  
  seqCodecDelete(codecp);
  seqFastqDelete(sqp);

  return errcode;
}
#endif

static int encodePAM(GUIDEWORD_T *wordF,
		     GUIDEWORD_T *maskF,
		     GUIDEWORD_T *wordR,
		     GUIDEWORD_T *maskR,
		     const char * const pamseq)
{
  int i, j, errcode = ERRCODE_SUCCESS;
  const int alphlen = (int) strlen(CODEC_ALPHABET);
  
  if (NULL == wordF || NULL == maskF ||
      NULL == wordR || NULL == maskR)
    return ERRCODE_NULLPTR;

  *wordF = *maskF = 0LL;
  *wordR = *maskR = 0LL;

  if (alphlen > UINT8_MAX)
    return ERRCODE_OVERFLOW;

  for (i=0; i<PAMSLEN; i++) {
    const int xf = NBITS_KEYCOD*(PAMSLEN-1-i);
    const int xr = NBITS_KEYCOD*(NBASES_WORD+i);
    const char pc = pamseq[i];
    if (pc == '\0')
      return ERRCODE_FAILURE;
    for (j = 0; j < alphlen && CODEC_ALPHABET[j] != pc; j++);
    if (j < SIZE_STANDARD_ALPHABET) { 
      *maskF |= ((GUIDEWORD_T) SEQCOD_STDNT_MASK)<<xf;
      *maskR |= ((GUIDEWORD_T) SEQCOD_STDNT_MASK)<<xr;
      *wordF |= ((GUIDEWORD_T) j)<<xf;
      *wordR |= ((GUIDEWORD_T) (SIZE_STANDARD_ALPHABET-1-j))<<xr;
    } else if (j >= alphlen) {
	return ERRCODE_FAILURE;
    }
  }
  if (pamseq[i] != '\0')
    errcode = ERRCODE_FAILURE;
  
  return errcode;
}

static int screenPAM(ErrMsg *errmsgp,
	      GuideIndex * const gxp,
	      const SeqSet * const ssp)
{
  int errcode = ERRCODE_SUCCESS;
  int nsu, minbasctr;
  char cod;
  unsigned char ssflg;
  //SEQNUM_t s, snum;
  SETSIZ_t siz;
  //SETSIZ_t const *soffsp;
  uint32_t const *basep;
  GUIDEPOS_T ctr, maxpos;
  GUIDEWORD_T word = 0LL, pamFmask = 0LL, pamRmask = 0LL;
  GUIDEWORD_T pamFkey=0LL, pamRkey = 0LL;

#ifdef DESIGUIDE_DEBUG
  SeqFastq *sqp = seqFastqCreate(0, SEQTYP_FASTA);
  SeqCodec *codecp = seqCodecCreate();

  if (NULL == sqp || NULL == codecp)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
#endif
  
  if (errcode != ERRCODE_SUCCESS)
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = encodePAM(&pamFkey, &pamFmask, &pamRkey, &pamRmask, 
			   gxp->pam)))
    ERRMSGNO(errmsgp, errcode);

  
  /* basep = (const uint32_t *) seqSetGetBaseData(&ssflg, &cod, &siz, ssp); */
  basep = (const uint32_t *) seqSetGetSeqDatByIndex(&ssflg, &cod, &siz, ssp);
  if (cod != SEQCOD_COMPRESSED || !(ssflg&SEQSET_COMPRESSED))
    ERRMSGNO(errmsgp, ERRCODE_SEQCODE);

  if (siz > GUIDEPOS_MAX)
    ERRMSGNO(errmsgp, ERRCODE_OVERFLOW);
  maxpos = (GUIDEPOS_T) siz;

  minbasctr = NBASES_WORD + PAMSLEN;
  for (nsu = MAXN_PER_UNIT-1, ctr = 0LL; ctr<maxpos; ctr++) {
    uint32_t bc = ((*basep)>>(nsu*NBITS_ALPHABET))&SEQCOD_ALPHA_MASK;
    if (--nsu < 0) {
      basep++;
      nsu = MAXN_PER_UNIT-1;
    }
    if (bc >= SIZE_STANDARD_ALPHABET) {
      minbasctr = NBASES_WORD + PAMSLEN;
    } else {
      word = (word<<NBITS_KEYCOD)|(bc&SEQCOD_STDNT_MASK);
      if (minbasctr > 0) {
	minbasctr--;
      } else {
	GUIDEPOS_T *posp;
	GUIDEWORD_T *guidep;
	int isFhit = 0LL == ((word^pamFkey)&pamFmask);
	int isRhit = 0LL == ((word^pamRkey)&pamRmask);
	if (isFhit) {
	  /* found a PAM sequence, add word and position to list */	
	  ARRNEXTP(guidep, gxp->guidr);
	  *guidep = (word>>(PAMSLEN*NBITS_KEYCOD))&GUIDE_MASK;
	  ARRNEXTP(posp, gxp->posr);
	  *posp = ctr - NBASES_WORD - PAMSLEN + 1;
	}
	if (isRhit) {
	  GUIDEWORD_T wordrc = word;
	  REVERSE_COMPLEMENT_WORD(wordrc, NBASES_WORD);
	  ARRNEXTP(guidep, gxp->guidr);
	  *guidep = (wordrc&GUIDE_MASK) | GUIDE_RCFLG;
	  ARRNEXTP(posp, gxp->posr);
	  *posp = ctr - NBASES_WORD + 1;
	}
#ifdef DESIGUIDE_DEBUG
	if ((isFhit) || (isRhit)) {
	  int i;
	  const char *seqstrp;
	  char cbuf[NBASES_WORD + PAMSLEN + 2];
	  SETSIZ_t os = ctr - NBASES_WORD - PAMSLEN + 1, oe = ctr;
	  /* check that the position and word are correct */
	  if ((errcode = seqSetFetchSegment(sqp, &os, &oe, ssp, codecp)))
	    ERRMSGNO(errmsgp, errcode);

	  seqstrp = seqFastqGetConstSequence(sqp, NULL, NULL);
	  decodeWord(cbuf, word, NBASES_WORD + PAMSLEN);
	  for (i=0; i<(NBASES_WORD + PAMSLEN) && cbuf[i] == seqstrp[i]; i++);
	  
	  /* printf("%s %llu %s %c%c\n",  */
	  /* 	 cbuf,  */
	  /* 	 (unsigned long long) ctr, */
	  /* 	 (i < NBASES_WORD + PAMSLEN)? "ERR": "OK", */
	  /* 	 (isFhit)?'F':' ', */
	  /* 	 (isRhit)? 'R':' '); */

	  if (i < NBASES_WORD + PAMSLEN)
	    ERRMSGNO(errmsgp, ERRCODE_SEQCODE);
	}
#endif
      }
    }
#ifdef DESIGUIDE_VERBOSE
    if ( 0 == (ctr & 0x7FFFFFF))
      printf("%5.2f %% processed ...\n", ((double) ctr)*100/siz);
#endif
  }
  printf("%5.2f %% processed ...\n", ((double) ctr)*100/siz);
  
#ifdef DESIGUIDE_DEBUG
  seqCodecDelete(codecp);
  seqFastqDelete(sqp);
#endif
  
  return errcode;
}

#ifdef DESIGUIDE_DEBUG
static int cmpGuideIndex(const GuideIndex * const ap, const GuideIndex * const bp)
{
  int i;
  size_t s, n = ARRLEN(ap->guidr);
  if (n != ARRLEN(bp->guidr) ||
      ARRLEN(ap->posr) != ARRLEN(bp->posr) ||
      n != ARRLEN(ap->posr))
    return ERRCODE_FAILURE;
  
  for (i=0; i<=PAMSLEN && ap->pam[i] == bp->pam[i]; i++);
  if (i <= PAMSLEN || ap->pam[PAMSLEN] != '\0' || bp->pam[PAMSLEN] != '\0')
    return ERRCODE_FAILURE;

  for (s = 0; s<n && ap->posr[s] == bp->posr[s]; s++);
  if (s < n)
    return ERRCODE_FAILURE;

  for (s = 0; s<n && ap->guidr[s] == bp->guidr[s]; s++);
  if (s < n)
    return ERRCODE_FAILURE;
    
  return ERRCODE_SUCCESS;
}
#endif

static int generateGuideIndex(ErrMsg *errmsgp,
			      const char * const fnamprfx_smaltidx,
			      const char * const pamseq,
			      const char * const fnam_idx)
{
  int errcode = ERRCODE_SUCCESS;
  SeqSet *ssp = NULL;
  GuideIndex guix;
#ifdef DESIGUIDE_DEBUG
  GuideIndex debugx;
#endif

  if ((errcode = initGuideIndex(&guix, MEMALLOC_BLKSZ)))
    ERRMSGNO(errmsgp, errcode);
  setGuideIndexPAM(&guix, pamseq);

  printf("loading genomic sequence ...\n");
  ssp = seqSetReadBinFil(&errcode, fnamprfx_smaltidx);
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  printf("recording PAM sequences ...\n");
  resetGuideIndex(&guix);
  if ((errcode = screenPAM(errmsgp, &guix, ssp)))
    ERRMSGNO(errmsgp, errcode);

  printf("writing index to files %s ...\n", fnam_idx);
  if ((errcode = writeGuideIndexToFile(&guix, 0, fnam_idx)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = writeGuideIndexToFile(&guix, 1, fnam_idx)))
    ERRMSGNO(errmsgp, errcode);

#ifdef DESIGUIDE_DEBUG 
  printf("checking index ...\n");
  if ((errcode = checkGuideIndex(&guix, ssp)))
     ERRMSGNO(errmsgp, errcode);

  printf("checking complete index ...\n");
  if ((errcode = checkGuideIndexComplete(&guix, ssp)))
     ERRMSGNO(errmsgp, errcode);

  printf("loading guide index ...\n");
  if ((errcode = initGuideIndex(&debugx, 0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if ((errcode = loadGuideIndexFromFile(&debugx, 0, fnam_idx)) ||
      (errcode = loadGuideIndexFromFile(&debugx, 1, fnam_idx)))
    ERRMSGNO(errmsgp, errcode);

  if (!(errcode = cmpGuideIndex(&guix, &debugx))) {
    printf("loaded guide index identical to orginal.\n");
  }
  freeGuideIndex(&debugx);
#endif

  seqSetDelete(ssp);
  freeGuideIndex(&guix);

  return errcode;
}

static int parseBEDfile(BED **bedpp, const char * const fnam_bed)
{
  int erc, errcode = ERRCODE_SUCCESS;
  size_t sl;
  char linbuf[LINBUF_MAX];
  BED *bp;
  FILE *fp;

  if (NULL == (fp = EFOPEN(fnam_bed, "r")))
    return ERRCODE_NOFILE;

  ARRNEXTP(bp, *bedpp);
  if (NULL == bp)
    return ERRCODE_NOMEM;

  while (ERRCODE_SUCCESS == errcode &&
	 NULL != fgets(linbuf, LINBUF_MAX, fp)) {
    int fldnum = 0, s = 0, e;
    int skipWhiteSpace = 1, skipline = 0;
    unsigned long uln;
    for (e = s; 
	 ERRCODE_SUCCESS == errcode &&
	   linbuf[e] != '\0'; 
	 e++) {
      if ((skipWhiteSpace)) {
	if (!isspace(linbuf[e])) {
	  s = e;
	  if (0 == fldnum && '#' == linbuf[e]) {
	    skipline = 1; /* comment/header line */
	    linbuf[e+1] = '\0';
	  }
	  skipWhiteSpace = 0;
	}
      } else if (isspace(linbuf[e])) {
	char *cp;
	fldnum++;
	linbuf[e] = '\0';
	switch (fldnum) {
	case 1: /* 1: chromosome name */
	  strncpy(bp->chrnam, linbuf + s, CHRNAM_MAX);
	  for (cp=bp->chrnam; *cp && !isspace(*cp); cp++);
	  *cp = '\0';
	  bp->genID[0] = '\0';
	  bp->uniprotID[0] = '\0';
	  bp->ensemblID[0] = '\0';
	  bp->annotation[0] = '\0';
	  break;
	case 2:
	  uln = strtoul(linbuf + s, (char **) NULL, 10);
	  if (ULONG_MAX == uln) {
	    errcode = ERRCODE_FILEFORM;
	    bp->start = 0;
	  } else {
	    bp->start = (SEQLEN_t) uln;
	  }
	  break;
	case 3:
	  uln = strtoul(linbuf + s, (char **) NULL, 10);
	  if (ULONG_MAX == uln) {
	    errcode = ERRCODE_FILEFORM;
	    bp->end = 0;
	  } else {
	    bp->end = (SEQLEN_t) uln;
	  }
	  break;
	case 4:	  
	  if ('+' == linbuf[s])
	    bp->strand = STRAND_FORWARD;
	  else if ('-' == linbuf[s])
	    bp->strand = STRAND_REVERSE_COMPLEMENT;
	  else
	    errcode = ERRCODE_FILEFORM;
	  strncpy(bp->annotation, linbuf + e + 1, FILELIN_MAX);
	  break;	  
	case 5:	  
	  strncpy(bp->genID, linbuf + s, GENID_MAX);
	  break;
	case 6:
	  break;
	case 7:
	  break;
	case 8:
	  strncpy(bp->uniprotID, linbuf + s, GENID_MAX);
	  break;
	case 9:
	  strncpy(bp->ensemblID, linbuf + s, GENID_MAX);
	  linbuf[e+1] = '\0'; /* ignore rest of fields */
	  break;
	default:
	  errcode = ERRCODE_FILEFORM;
	}
	skipWhiteSpace = 1;
      }
    }
    if ((fldnum >= 4) &&
	((isspace(linbuf[e-1])) || linbuf[e-1] == '\0')) {
      ARRNEXTP(bp, *bedpp);
    } else if (ERRCODE_SUCCESS == errcode && !skipline) {
      errcode = ERRCODE_FILEFORM;
    }
  }
  if (strlen(bp->uniprotID) < 1)
    strcpy(bp->uniprotID, NOSTR);
  sl = strlen(bp->ensemblID);
  if (sl > 0 && isspace(bp->ensemblID[sl-1])) {
    sl--;
    bp->ensemblID[sl] = '\0';
  }
  if (sl < 1)
    strcpy(bp->ensemblID, NOSTR);
  erc = EFCLOSE(fp);
  ARRLEN(*bedpp)--;

  return ((errcode))? errcode: erc;
}

static int addSeqSetInfoToBED(BED * const bedr, const SeqSet * const ssp)
{
  size_t b, nbed = ARRLEN(bedr);
  const SETSIZ_t *soffsp;
  SEQNUM_t s;
  const SEQNUM_t nseq = seqSetGetOffsets(ssp, &soffsp);
  const char **namr = NULL;
  
  if (NULL == ECALLOCP(nseq, namr))
    return ERRCODE_NOMEM;

  for (s=0; s<nseq; s++) {
    seqSetGetSeqDatByIndex(NULL, namr + s, s, ssp);
  }
  
  for (b=0; b<nbed; b++) {
    char *cnam = bedr[b].chrnam;
    for (s=0; s<nseq; s++) {
      char *np = cnam;
      const char *cp = namr[s];
      for (; *cp == *np && *np != '\0'; cp++, np++);
      if (*np == '\0' && (*cp == '\0' || isspace(*cp)))
	break;
    }
    if (s >= nseq)
      return ERRCODE_FAILURE;
    bedr[b].seqidx = s;
    bedr[b].os = soffsp[s] + bedr[b].start - 1;
    bedr[b].oe = soffsp[s] + bedr[b].end - 1;
  }
  
  free(namr);
  return ERRCODE_SUCCESS;
}

static int writeReport(FILE *fp, int with_offtargets, 
		       const GuideIndex * const gxp, uint32_t idx, 
		       const MATCH * const matchr, int n_mismatch_min,
		       const BED * const bedp,
		       SeqFastq * const sqp,
		       const SeqSet * const ssp, const SeqCodec * const codecp)
{
  int isRC, n_char, n_offtarget, n_score_offtarget_min, score_offtarget_min, errcode;
  uint32_t i, nmm_report, nm = ARRLEN(matchr), ng = ARRLEN(gxp->guidr);
  char cbuf[NBASES_WORD + 1];
  char pambuf[PAMSLEN + 1];
  char refbuf[NBASES_WORD + PAMSLEN + 1];
  const char *seqstrp;
  GUIDEWORD_T word;
  GUIDEPOS_T pos;
  SEQLEN_t slen;
  SETSIZ_t chrom_offs, set_start, set_end;
  const MATCH *matchptrs[OFFTARGETS_MAXNUM_REPORT];
  
  if (ARRLEN(gxp->guidr) > UINT32_MAX || 
      ARRLEN(matchr) > UINT32_MAX) 
    return ERRCODE_OVERFLOW;

  if (idx > ng)
    return ERRCODE_ARGRANGE;

  word = gxp->guidr[idx];
  pos = gxp->posr[idx];
  isRC = (word & GUIDE_RCFLG) != 0LL;

  if ((isRC)) {
    set_start = pos - PAMSLEN;
  } else {
    set_start = pos;
  }
  set_end = set_start + NBASES_WORD + PAMSLEN - 1;
  chrom_offs = set_start - bedp->os + bedp->start;

  if ((errcode = seqSetFetchSegment(sqp, &set_start, &set_end, 
				    ssp, codecp)))
    return errcode;
 
  seqstrp = seqFastqGetConstSequence(sqp, &slen, NULL);
  if (slen != NBASES_WORD + PAMSLEN)
    return ERRCODE_FAILURE;

  strncpy(refbuf, seqstrp, NBASES_WORD + PAMSLEN);
  refbuf[NBASES_WORD + PAMSLEN] = '\0';

  if ((isRC) && (errcode = seqFastqReverse(sqp, codecp)))
    return errcode;
  seqstrp = seqFastqGetConstSequence(sqp, &slen, NULL);
  strncpy(pambuf, seqstrp + NBASES_WORD, PAMSLEN);
  pambuf[PAMSLEN] = '\0';

  /* establish number of off-target hits with the same number of n_mismatch_min
     mismatches */
  n_offtarget = 0;
  score_offtarget_min = NBASES_WORD*NBASES_WORD;
  n_score_offtarget_min = 0;
  for (i = 0; i<nm; i++) {
    if (matchr[i].nmm <= n_mismatch_min) {
      if (n_offtarget < OFFTARGETS_MAXNUM_REPORT) {
	const MATCH * const mp = matchr + i;
	matchptrs[n_offtarget] = mp;
	const int sc = scoreMatch(NULL, mp);
	if (sc <= score_offtarget_min) {
	  if (sc == score_offtarget_min) {
	    n_score_offtarget_min++;
	  } else {
	    score_offtarget_min = sc;
	    n_score_offtarget_min = 1;
	  }
	}
      }
      n_offtarget++;
    }
  }
  nmm_report = (n_offtarget <= OFFTARGETS_MAXNUM_REPORT)?
    n_offtarget: OFFTARGETS_MAXNUM_REPORT;
   
 if ((with_offtargets))
   fprintf(fp, "TARGET ");
 /* [1] location string of target segment (e.g.exon), example: 2:164769979-164770080;+
  * [2] guide seqeuence
  * [3] PAM sequence along reference (forward) strand
  * [4] (reverse complement of) matched reference segment
  * [5] 1-based start of match (including PAM) on reference
  * [6] 1-based end of match (including PAM) on reference
  * [7] matched strand [+-]
  * [8] smallest number of mismatches for guide among off-targets
  * [9] number of off-targets with this number of mismatches
  * [10] smallest score among off-targets
  * [11] number of off-targets with that minimum score
  * [12 ...] bed file annotation fields, i.e. HGNC gene symbol, ENSEMBL gene ID, ....
  */
  /* n_char = fprintf(fp, "%s\t%s\t%s\t%s:%u-%u;%c\t%s\t%s\t%s\t%llu\t%llu\t%c\t%i\t%i", */
  /* 		   bedp->uniprotID, bedp->genID, bedp->ensemblID,  */
  /* 		   bedp->chrnam, bedp->start, bedp->end,  */
  /* 		   (bedp->strand < 0)? '-':'+', */
  /* 		   decodeWord(cbuf, word, NBASES_WORD), */
  /* 		   pambuf, refbuf,  */
  /* 		   (ULL_T) chrom_offs + 1, (ULL_T) chrom_offs + NBASES_WORD +PAMSLEN, */
  /* 		   (isRC)? '-':'+', */
  /* 		   n_mismatch_min, n_offtarget); */
  n_char = fprintf(fp, "%s:%u-%u;%c\t%s\t%s\t%s\t%llu\t%llu\t%c\t%i\t%i\t%i\t%i\t%s",
		   bedp->chrnam, bedp->start, bedp->end, 
		   (bedp->strand < 0)? '-':'+',
		   decodeWord(cbuf, word, NBASES_WORD),
		   pambuf, refbuf, 
		   (ULL_T) chrom_offs + 1, (ULL_T) chrom_offs + NBASES_WORD + PAMSLEN,
		   (isRC)? '-':'+',
		   n_mismatch_min, n_offtarget,
		   score_offtarget_min, n_score_offtarget_min,
		   bedp->annotation
		   );

 
  if ((with_offtargets)) {
    for (i = 0; i < nmm_report && n_char > 0; i++) {
      const char *chrnam;
      int c;
      char chrnambuf[CHRNAM_MAX+1];
      char mmstr[NBASES_WORD + 1];
      uint32_t idx = matchptrs[i]->idx;
      int offt_score = scoreMatch(mmstr, matchptrs[i]);
      GUIDEWORD_T offt_word = gxp->guidr[idx];
      GUIDEPOS_T offt_pos = gxp->posr[idx];
      int offt_isRC = (offt_word & GUIDE_RCFLG) != 0LL;
      SEQNUM_t seqidx;
      SETSIZ_t seqoffs;
      SEQLEN_t so;
      SETSIZ_t offs = (isRC)? offt_pos - PAMSLEN: offt_pos;
      
      if ((errcode = seqSetGetIndexAndOffset(&seqidx, &seqoffs, offs, ssp)))
	return errcode;

      seqSetGetSeqDatByIndex(NULL, &chrnam, seqidx, ssp);
      for (c=0; c<CHRNAM_MAX && chrnam[c] && !isspace(chrnam[c]); c++)
	chrnambuf[c] = chrnam[c];    
      chrnambuf[c] = '\0';

      so = offs - seqoffs;
      n_char = fprintf(fp, "\nOFFTARGET[%i]\t%s\t%s:%u-%u:%c\t%i\t%s", i,
		       decodeWord(cbuf, offt_word, NBASES_WORD),
		       chrnambuf, so+1, so + NBASES_WORD,
		       (offt_isRC)? '-':'+',
		       offt_score, mmstr);
    }
  }
  fprintf(fp, "\n");
  return (n_char > 0)? ERRCODE_SUCCESS: ERRCODE_WRITEERR;
}

static int reportMatchesForGuide(FILE *fp, 
				 SeqFastq * const sqbufp,
				 const GuideIndex * const gxp,
				 const MATCH * const matchr,
				 int nmm_cutoff,
				 int nmm_min,
				 int margins[2],
				 const SeqFastq * const sqp,
				 const SeqSet * const ssp,
				 const SeqCodec * const codecp)
{
  int errcode = ERRCODE_SUCCESS;
  const uint32_t nmatch = ARRLEN(matchr);
  char nambuf[NAMBUFSZ];
  char cbuf[NBASES_WORD + 1];
  char gbuf[NBASES_WORD + 1];
  char mbuf[NBASES_WORD + 1];
  uint32_t i, linctr;
  int nmm, nmm_check;
  const char * const guidnam = seqFastqGetSeqName(sqp);
  if (ARRLEN(gxp->guidr) > UINT32_MAX || 
      ARRLEN(matchr) > UINT32_MAX) 
    return ERRCODE_OVERFLOW;


  if (nmatch < 1) 
    return ERRCODE_SUCCESS;
 
  if (nmm_min > nmm_cutoff) 
    nmm_cutoff = nmm_min;
  
  linctr = 1;
  for (nmm = nmm_min; nmm <= nmm_cutoff; nmm++) {
    for (i=0; i<nmatch; i++) {
      if  (matchr[i].nmm == nmm) {
	int mmscore;
	uint32_t idx = matchr[i].idx;
	GUIDEWORD_T word = gxp->guidr[idx];
	GUIDEPOS_T pos = gxp->posr[idx];
	/* leftmost position of the word */
	SEQNUM_t seqidx;
	SEQLEN_t so, rs, m1, m2;
	SEQLEN_t sl = NBASES_WORD;
	SEQLEN_t m5P = (margins[0] > 0)? (SEQLEN_t) margins[0]: 0;
	SEQLEN_t m3P = (margins[1] > 0)? (SEQLEN_t) margins[1]: 0;
	SETSIZ_t seqoffs;
	char strand = (word & GUIDE_RCFLG)? '-': '+';
	const char *seqnam;
	char *seqstrp;
	if ((errcode = seqSetGetIndexAndOffset(&seqidx, NULL, pos, ssp)))
	  return errcode;
	seqSetGetSeqDatByIndex(&seqoffs, &seqnam, seqidx, ssp);
	so = pos - seqoffs;
	getShortSequenceName(nambuf, NAMBUFSZ-1, seqnam);
	decodeWord(cbuf, word, NBASES_WORD);
	decodeWord(gbuf, matchr[i].gw, NBASES_WORD);
	seqFastqBlank(sqbufp);

	if ((word & GUIDE_RCFLG)) {
	  m1 = m3P;
	  m2 = m5P;
	} else {
	  m1 = m5P;
	  m2 = m3P;	    
	}
	if (so > m1) {
	  sl += m1 + m2;
	  rs = so - m1;
	} else {
	  sl += m2 + so;
	  rs = 0;
	}
	errcode = seqSetFetchSegmentBySequence(sqbufp, seqidx, rs, sl,
					       ssp, codecp);
	if ((errcode)) 
	  return errcode;

	if ((word & GUIDE_RCFLG) && (errcode = seqFastqReverse(sqbufp, codecp)))
	  return errcode;
	
	seqstrp = seqFastqGetSequence(sqbufp, NULL, NULL);
	nmm_check = generateMisMatchString(mbuf, &mmscore, gbuf, seqstrp); 
	/* fprintf(fp, "%s [%i] %u %llu %s %lu %c %i %s\n", guidnam,
	 *    idx, (unsigned int) seqidx, (unsigned long long) pos,
	 *    nambuf, (unsigned long) so, strand, nmm, cbuf);
	 */
	fprintf(fp, "%s\t%s\t[%i]\tchr%s:%lu-%lu:%c\t%s\t%i\t%s\t%s\t%i", 
		guidnam, gbuf, linctr, 
		nambuf, (unsigned long) so + 1, 
		(unsigned long) so + NBASES_WORD,
		strand, cbuf, nmm, seqstrp, mbuf, mmscore);
	convertSeq2RY(seqstrp);
	fprintf(fp, "\t%s\t%s\n", seqstrp, (nmm_check == nmm)? "OK": "ERR");
	linctr += 1;
      }
    }
  }
  
  return errcode;
}

static int findWordInGuides(const GUIDEWORD_T testword,
			    const GUIDEWORD_T * const guidr,
			    int n_mismatch_cutoff,
			    int *n_mismatch_min, MATCH **matchrp)
{
  int nmm_min = (int) (sizeof(GUIDEWORD_T)*8/NBITS_KEYCOD); /* initialised to maximum */
  int nmm_stop = nmm_min;
  const uint32_t ng = (uint32_t) ARRLEN(guidr);
  uint32_t i;

  for (i=0; i<ng; i++) {
    MATCH *mp;
    int nmm;
    GUIDEWORD_T const gw = (guidr[i]^testword) & GUIDE_MASK;
    GUIDEWORD_T df;
    for (df = gw, nmm=0; df != 0LL && nmm <= nmm_stop; df >>= NBITS_KEYCOD) {
      if ((df&SEQCOD_STDNT_MASK))
	nmm++;
    }
    if (nmm > nmm_stop)
      continue;

    ARRNEXTP(mp, *matchrp);
    if (NULL == mp) 
      return ERRCODE_NOMEM;
    mp->idx = i;
    mp->df = gw;
    mp->gw = testword;
    mp->nmm = nmm;
    if (nmm < nmm_min) {
      nmm_min = nmm;
      if (nmm_min < n_mismatch_cutoff)
	nmm_stop = n_mismatch_cutoff;
      else
	nmm_stop = nmm_min;
    }
  }

  if (n_mismatch_min != NULL) *n_mismatch_min = nmm_min;
  
  return ERRCODE_SUCCESS;
}

static int scanGuides(const GUIDEWORD_T * const guidr, size_t x, 
		      int *n_mismatch_min, MATCH **matchrp)
{
  int nmm_min = (int) sizeof(GUIDEWORD_T)*4;
  const GUIDEWORD_T guidemask = ~(~((GUIDEWORD_T) 0) << NBITS_KEYCOD*SEEDCORE);
  uint32_t i;
  uint32_t ng = (uint32_t) ARRLEN(guidr);
  const GUIDEWORD_T testword = guidr[x];

  if (n_mismatch_min != NULL) *n_mismatch_min = 0;
  ARRLEN(*matchrp) = 0LL;

  if (ARRLEN(guidr) > UINT32_MAX)
    return ERRCODE_OVERFLOW;

  if (x >= ng)
    return ERRCODE_ARGRANGE;

  for (i=0; i<ng; i++) {
    MATCH *mp;
    
    GUIDEWORD_T df = (guidr[i]^testword) & guidemask;
    int nmm = 0;
    
    while (df != 0LL && nmm <= nmm_min) {
      if ((df&SEQCOD_STDNT_MASK))
	nmm++;     
      df >>= NBITS_KEYCOD;
    }
    if (nmm > nmm_min || x == i)
      continue;
    ARRNEXTP(mp, *matchrp);
    if (NULL == mp) 
      return ERRCODE_NOMEM;
    mp->idx = i;
    mp->df = (guidr[i]^testword) & GUIDE_MASK;
    mp->gw = guidr[i];
    mp->nmm = nmm_min = nmm;
    if (0 == nmm_min) 
      break;
  }

  if (n_mismatch_min != NULL) *n_mismatch_min = nmm_min;

  return ERRCODE_SUCCESS;
}

static int 
selectGuidesInSetOfBedSegments
( ErrMsg * const errmsgp,
  const BED * const bedr,
  const GuideIndex * const guixp,
  const int edit_distance_min,
  FILE * const oufh_guides,
  FILE * const oufh_filtered,
  FILE * const oufh_offtarget,
  MATCH ** const matchrp,
  SeqFastq * const sqbufp,
  const SeqSet * const ssp,
  const SeqCodec * const codecp)
{
  int isOK, errcode = ERRCODE_SUCCESS;
  const uint32_t margin = NBASES_WORD + PAMSLEN;
  uint32_t g = 0, b = 0;
  const uint32_t nbed = ARRLEN(bedr);
  const uint32_t nguides = ARRLEN(guixp->guidr);

  isOK = g<nguides && b < nbed;
  while (isOK) {
    if (guixp->posr[g] < bedr[b].os) {
      g++;
      //printf("guide %i (%6.3f %%) ...\n", g, ((double) (g*100))/nguides);
      isOK = g<nguides;
    } else if (guixp->posr[g] > bedr[b].oe + margin) {
      b++;
      isOK = b < nbed;
      printf("%u (%6.3f %%) of BED segments completed ...\n", 
	     b, ((double) (b*100))/nbed); 
      if ((isOK) && g>0 && guixp->posr[g-1] > bedr[b].os) { 
	/* successive bed segments may overlap -> scroll guides back 
	 * to beginning of next bed segment if necessary */
	printf("BED[%u]: scrolling back from %u ...\n", b, g);
	for (; g>0 && guixp->posr[g-1] > bedr[b].os; g--);
	printf("         to %u.\n", g);
      }
      fflush(oufh_guides);
      fflush(oufh_filtered);
      fflush(oufh_offtarget);
    } else {
      int nmm_min;
      if ((errcode = scanGuides(guixp->guidr, g, &nmm_min, matchrp))) {
	ERRMSGNO(errmsgp, errcode);
      	isOK = 0;
      } else {
      	FILE *oufp = (nmm_min >= edit_distance_min)?
      	  oufh_guides: oufh_filtered;
      	writeReport(oufp, 0, guixp, g, *matchrp, nmm_min, bedr+b,
      		    sqbufp, ssp, codecp);
      	if (nmm_min >= edit_distance_min) {
      	  writeReport(oufh_offtarget, 1, guixp, g, *matchrp, nmm_min, bedr+b,
      		      sqbufp, ssp, codecp);
      	}
	g++;
	isOK = g<nguides;
	printf("guide %i (%6.3f %%) ...\n", g, ((double) (g*100))/nguides);
      }
    }
  }
  return errcode;
}

/* static int  */
/* selectGuidesFromSingleBedSegment */
/* ( ErrMsg * const errmsgp, */
/*   const BED * const bedp, /\**< BED segment *\/ */
/*   const GuideIndex * const guixp, */
/*   uint32_t *gx,  */
/*   /\**< index for starting scrolling in GuideIndex, gets updated *\/ */
/*   const int edit_distance_min, */
/*   FILE * const oufh_guides, */
/*   FILE * const oufh_filtered, */
/*   FILE * const oufh_offtarget, */
/*   MATCH ** const matchrp, */
/*   SeqFastq * const sqbufp, */
/*   const SeqSet * const ssp, */
/*   const SeqCodec * const codecp) */
/* { */
/*   int errcode = ERRCODE_SUCCESS; */
/*   const uint32_t nguides = ARRLEN(guixp->guidr); */
/*   const uint32_t margin = NBASES_WORD + PAMSLEN; */
/*   uint32_t x = (NULL == gx)? 0: *gx; */
  
/*   /\* scroll to first base upstream *\/ */
/*   if (bedp->strand < 0) { */
/*     for (; guixp->posr[x] > bedp->oe && x>0; x--); */
/*     for (; guixp->posr[x] < bedp->oe && x < nguides; x++); */
/*   } else { */
/*     for (; guixp->posr[x] + margin < bedp->os && x < nguides; x++); */
/*     for (; guixp->posr[x] + margin > bedp->os && x > 0; x--); */
/*   } */
  
/*   while (0) { */
    
/*   } */
/* } */

static int selectGuidesFromBED(ErrMsg *errmsgp,
			       const int edit_distance_min,
			       const char * const fnamprfx_smaltidx,
			       const char * const fnamprfx_guidx,
			       const char * const fnam_bed,
			       const char * const fnamprfx_sel)
{
  int errcode = ERRCODE_SUCCESS;
  char fnambuf[FILENAME_MAX];
  uint32_t nbed;
  BED *bedr;
  MATCH *matchr;
  GuideIndex guix;
  SeqFastq *sqp = seqFastqCreate(0, SEQTYP_FASTA);
  SeqCodec *codecp = seqCodecCreate();
  SeqSet *ssp;
  FILE *oufh_guides, *oufh_filtered, *oufh_offtarget;

  oufh_guides = 
    EFOPEN(strcat(strcpy(fnambuf, fnamprfx_sel), "_guides.csv"), "w");
		       
  oufh_filtered = 
    EFOPEN(strcat(strcpy(fnambuf, fnamprfx_sel), "_filtered.csv"), "w");

  oufh_offtarget = 
    EFOPEN(strcat(strcpy(fnambuf, fnamprfx_sel), "_offtargets.csv"), "w");

  if (NULL == oufh_guides ||
      NULL == oufh_filtered ||
      NULL == oufh_offtarget)
    return ERRCODE_FILEIO;

  printf("loading genomic sequence ...\n");
  ssp = seqSetReadBinFil(&errcode, fnamprfx_smaltidx);
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  printf("loading guide index ...\n");
  if ((errcode = initGuideIndex(&guix, 0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if ((errcode = loadGuideIndexFromFile(&guix, 0, fnamprfx_guidx)) ||
      (errcode = loadGuideIndexFromFile(&guix, 1, fnamprfx_guidx)))
    ERRMSGNO(errmsgp, errcode);

#ifdef DESIGUIDE_DEBUG
  printf("checking guide index ...\n");
  if ((errcode = checkGuideIndex(&guix, ssp)))
    ERRMSGNO(errmsgp, errcode);
  else
    printf("guide index OK.\n");
#endif

  printf("loading BED file ...\n");
  if (NULL == ARRCREATE(bedr, 0))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((errcode = parseBEDfile(&bedr, fnam_bed)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = addSeqSetInfoToBED(bedr, ssp)))
     ERRMSGNO(errmsgp, errcode);

  nbed = ARRLEN(bedr);  
  printf("sorting %u BED segments ...\n", nbed);
  
  qsort(bedr, (size_t) nbed, sizeof(BED), cmpBED);

  if (ARRLEN(bedr) > UINT32_MAX ||
      ARRLEN(guix.guidr) > UINT32_MAX ||
      ARRLEN(guix.guidr) != ARRLEN(guix.posr))
    ERRMSGNO(errmsgp, ERRCODE_OVERFLOW);

  if (NULL == ARRCREATE(matchr, 0))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  printf("select guides in BED segments ...\n");
  errcode = selectGuidesInSetOfBedSegments(errmsgp, 
					   bedr, &guix, 
					   edit_distance_min,
					   oufh_guides, oufh_filtered,
					   oufh_offtarget,
					   &matchr,
					   sqp, ssp, codecp);
  ARRDELETE(matchr);
  ARRDELETE(bedr);

  seqFastqDelete(sqp);
  seqCodecDelete(codecp);

  freeGuideIndex(&guix);
  seqSetDelete(ssp);

  EFCLOSE(oufh_offtarget);
  EFCLOSE(oufh_filtered);
  EFCLOSE(oufh_guides);

  return errcode;
}

static int 
findGuideMatchesInIndex
(ErrMsg *errmsgp,
 FILE *oufp,
 SeqFastq * const sqbufp,
 const SeqFastq * const sqp,
 int nmm_cutoff,
 int margins[2],
 const GuideIndex * const guixp,
 const SeqSet * const ssp,
 const SeqCodec * const codecp) 
{
  int errcode = ERRCODE_SUCCESS;
  GUIDEWORD_T guide = 0LL;
  int nmm_min = 0;
  MATCH *matchr;
#ifdef DESIGUIDE_DEBUG
  char cbuf[NBASES_WORD + 1];
#endif
  if ((errcode = encodeWord(&guide, sqp)))
      ERRMSGNO(errmsgp, errcode);

#ifdef DESIGUIDE_DEBUG
  decodeWord(cbuf, guide, NBASES_WORD);
  fprintf(oufp, "DECODED %s\n", cbuf);
#endif

  if (NULL == ARRCREATE(matchr, 0))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((errcode = findWordInGuides(guide, guixp->guidr, nmm_cutoff,
				  &nmm_min, &matchr)))
     ERRMSGNO(errmsgp, errcode);

  if ((errcode = reportMatchesForGuide(oufp, sqbufp, guixp, matchr,
				       nmm_cutoff, nmm_min,
				       margins,
				       sqp, ssp, codecp)))
      ERRMSGNO(errmsgp, errcode);

  ARRDELETE(matchr);

  return errcode;
}

static int 
checkGuidesForOfftargets
(ErrMsg *errmsgp, 
 const char * const fnam_output,
 const int nmm_cutoff,
 int margins[],
 const char * const fnamprfx_smaltidx,
 const char * const fnamprfx_guidx, 
 const char * const fnam_guides)
{
  int errcode = ERRCODE_SUCCESS;
  SEQNUM_t sctr;
  FILE *oufp = NULL;
  SeqFastq *sqp, *sqbufp;
  SeqCodec *codecp;
  SeqIO *sfp;
  SeqSet *ssp;
  GuideIndex guix;
  
  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sqp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(sqbufp = seqFastqCreate(0, SEQTYP_FASTA)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  

  printf("loading genomic sequence ...\n");
  ssp = seqSetReadBinFil(&errcode, fnamprfx_smaltidx);
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  printf("loading guide index ...\n");
  if ((errcode = initGuideIndex(&guix, 0)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if ((errcode = loadGuideIndexFromFile(&guix, 0, fnamprfx_guidx)) ||
      (errcode = loadGuideIndexFromFile(&guix, 1, fnamprfx_guidx)))
    ERRMSGNO(errmsgp, errcode);

#ifdef DESIGUIDE_DEBUG
  printf("checking guide index ...\n");
  if ((errcode = checkGuideIndex(&guix, ssp)))
    ERRMSGNO(errmsgp, errcode);
  else
    printf("guide index OK.\n");
#endif

  if (NULL == (sfp = seqIOopen(&errcode, fnam_guides, SEQIO_READ, 0)))
    ERRMSGNO(errmsgp, errcode);

  if (!(oufp = EFOPEN(fnam_output, "wb")))
    return ERRCODE_NOFILE;
 
  fprintf(oufp, "guide_name\tguide_sequence\t[running_index]\tchromosome:start-end:strand"\
	  "\ttarget_sequence\tn_mismatches\tref_sequence\tRY_seq\n");
  sctr = 0;
  while(!seqIOstatus(sfp)) {
    printf("Reading query sequence %li ...\n", sctr++);
    if ((errcode = seqFastqRead(sqp, sfp)))
      ERRMSGNO(errmsgp, errcode);

#ifdef DESIGUIDE_DEBUG
    fprintf(oufp, "%s %s\n", 
	    seqFastqGetSeqName(sqp), 
	    seqFastqGetConstSequence(sqp, NULL, NULL));
#endif
    if ((errcode = seqFastqEncode(sqp, codecp)))
       ERRMSGNO(errmsgp, errcode);

    if ((errcode = seqFastqCompress(sqp)))
       ERRMSGNO(errmsgp, errcode);

    if ((errcode = findGuideMatchesInIndex(errmsgp, oufp,
					   sqbufp, sqp, 
					   nmm_cutoff, margins,
					   &guix, ssp, codecp)))
      ERRMSGNO(errmsgp, errcode);
    fflush(oufp);
  }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sfp));

  EFCLOSE(oufp);
  seqIOclose(sfp);
  freeGuideIndex(&guix);
  seqSetDelete(ssp);
  seqFastqDelete(sqbufp);
  seqFastqDelete(sqp);
  seqCodecDelete(codecp);

  return errcode;
}

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  //time_t time_start, time_stop;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc == 5 && 0 == strcmp(argv[1], "index")) {
    char *fnamprfx_smaltidx = argv[2];
    char *pamseq = argv[3];
    char *fnam_idx = argv[4];
    generateGuideIndex(errmsgp, fnamprfx_smaltidx, pamseq, fnam_idx);
  } else if (argc == 6 && 0 == strcmp(argv[1], "select")) {
    char *fnamprfx_smaltidx = argv[2];
    char *fnamprfx_guidx = argv[3];
    char *fnam_bed = argv[4];
    char *fnamprfx_sel = argv[5];

    errcode = selectGuidesFromBED(errmsgp, 
				  MINIMUM_EDIT_DISTANCE,
				  fnamprfx_smaltidx,
				  fnamprfx_guidx, fnam_bed, fnamprfx_sel);
  } else if ((argc == 7 || argc == 8) && 0 == strcmp(argv[1], "check")) {
    char *fnamprfx_smaltidx = argv[2];
    char *fnamprfx_guidx = argv[3];
    int nmm_cutoff = atoi(argv[4]);
    char *fnam_guides_fasta = argv[5];
    char *fnam_output = argv[6];
    int margins[2] = {0,0};

    if (argc > 7) {
      char *as = argv[7];
      int i, sl = (int) strlen(as);
      for (i=0; i<sl && as[i] != ','; i++);
      as[i] = '\0';
      margins[0] = atoi(as);
      if (i < sl) {
	margins[1] = atoi(as+i+1);
      }
    }
    checkGuidesForOfftargets(errmsgp, fnam_output, nmm_cutoff,
			     margins,
			     fnamprfx_smaltidx, fnamprfx_guidx, 
			     fnam_guides_fasta);  
  } else {
    fprintf(stderr, 
	    "usage: %s index"\
	    " <SMALT index prefix> <PAM sequence> <guide index prefix>\n" \
	    "       %s select"\
	    " <SMALT index prefix> <guide index prefix> <BED file>"\
	    " <guide index prefix for selected guides (output)>\n"\
	    "       %s check"\
	    " <SMALT index prefix> <guide index prefix>"\
	    " <max num mismatches> <guide sequences (FASTA)>"\
	    " <output file> [<5'-margin>,<3'-margin>]\n",
	    argv[0], argv[0], argv[0]);

    return ERRCODE_FAILURE;
  }

  printf("finished.\n");
  
  ERRMSG_END(errmsgp);
  return errcode;
}
