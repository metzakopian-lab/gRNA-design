#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdio>
#include <cstring>

#include "Model.h"
#include "Serialization.h"

#define BUFFER_SIZE 120000

#define N_STATE 1
#define EOF_STATE 0
#define SEQ_STATE 1 
#define CHROMOSOME_STATE 2

GuideModel gRNAs;

int readFastaLine(FILE *fp, char* buffer, const size_t buf_size, size_t* length)
{
  char c;
  size_t i;
  for ( i = 0; i < buf_size && (c = fgetc(fp)) != '>' && c != EOF; ++i)
  {
    if(c == '\n')
    {
      continue;
    }
    buffer[i] = c;
  }

  buffer[i] = '\0';
  
  *length = i;

  if (c == EOF)
  {
    return EOF_STATE;
  }
  else if (c == '>')
  {
    return CHROMOSOME_STATE;
  }
  else
  {
    return SEQ_STATE;
  }
}

void extractgRNA(const char* buf, size_t new_bytes, int nt, const std::string& pam, const std::string& chr_name, size_t& genome_offset_counter, size_t& guide_uid)
{
  
  char gRNA[256];
  unsigned offset = nt + pam.size();
  const char* end = buf + new_bytes + offset;
  unsigned char pam_state = 0;
  
  for (const char* b = buf + nt; b != end; ++b, ++genome_offset_counter)
  {
    if (pam[pam_state] == 'N')
    {
      pam_state++;
    }
    else if (pam[pam_state] == *b)
    {
      pam_state++;
    }
    else
    {
      pam_state = 0;
    }
    if(pam_state == pam.size())
    {
      memcpy(gRNA,b - offset, nt);
      gRNA[nt] = '\0';
      bool all_nts = true;
      for(int i = 0; i < nt; i++)
      {
        const char& c = gRNA[i];
        if( not ((c == 'A') || (c == 'T') || (c == 'G') || (c == 'G')))
        {
          all_nts = false;
          break;
        }
      }
      if (all_nts)
      {
        GuideMeta g;
        g.start_pos = genome_offset_counter - pam.size() - nt;
        g.end_pos = g.start_pos + nt;
        g.chromosome = chr_name;
        g.id = guide_uid++;
        auto it = gRNAs.insert(std::make_pair(gRNA, std::vector<GuideMeta>()));
        it.first->second.push_back(g);
      }
    }
  }
}

int worker(FILE* fp, unsigned char nt, const std::string pam)
{
  int write_buf_pos_index = 0;
  const size_t buf_size = 100000000;
  char* buffer = new char[buf_size]; // 10M chunks

  size_t buf_read = 0;
  char* buf = buffer;


  size_t offset = nt + pam.size();
  size_t pos_counter = 0;
  size_t guide_id = 0;
  std::string chr;

  while (1)
  {
    int state = readFastaLine(fp, buf, buf_size - offset, &buf_read);
    std::string chr_name;
    if (state == CHROMOSOME_STATE)
    {
      // read the new chromosome
      char chromname[1024];
      fgets(chromname, 1024, fp);


      chr = chromname;
      
      std::size_t pos = chr.find_first_of(" \t");
      if(pos != std::string::npos)
      {
        chr_name = chr.substr(0,pos);
      }
      
      std::cout << " Processing chromosome: " << chromname << " as " << chr_name << std::endl;
      buf = buffer;
      pos_counter = 0;
      continue;

    }
    extractgRNA(buffer, buf_read, nt, pam, chr_name, pos_counter, guide_id);
    memcpy(buffer, buf + buf_read - offset, offset);
    buf = buffer + offset;
    
    if (state == EOF_STATE)
    {
      // Let's get the fuck out of here
      break;
    }
  }
  
  return 1;
}



int main(int argc, char* argv[])
{
  
  if(argc != 3)
  {
    std::cerr << "Invalid number of arguments" << std::endl;
    std::cout << "Usage: " << argv[0] << "<fasta file> <binary-archive>" << std::endl;
    return 1;
  }

  std::string file = argv[1];
  std::string gRNA_file = argv[2];
  
  std::cout << "Input file " << file << std::endl;
  FILE *fp = fopen(file.c_str(), "r");
  
  if (fp)
  {
    worker(fp,19,"NGG");
  }
  else
  {
    std::cerr << file << " not found" << std::endl;
    return 1;
  }
  size_t t = 0, max = 0;
  
  modelSerialize(gRNAs,gRNA_file);
 
  return 0;
}


