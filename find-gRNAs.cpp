#include <iostream>
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include "Tree.h"

#define BUFFER_SIZE 120000
#define CHROMOSOME_STATE 2
#define N_STATE 1

std::string current_chromosome;
int read_pos = 0;
int write_pos = 0;
int readFastaLine(FILE *fp, char* buffer, size_t buf_size, size_t* length)
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
    return 0;
  }
  else if (c == '>')
  {
    return 2;
  }
  else
  {
    return 1;
  }

}

typedef struct{
  size_t start_pos;
  size_t end_pos;
  size_t id;
  std::string chromosome;


}gRNAmeta;

std::map<std::string, gRNAmeta> gRNAs;

void extractgRNA(const char* buf, size_t new_bytes, int nt, const std::string& pam, const std::string& chr_name, size_t& genome_offset_counter)
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
      memcpy(gRNA, b - offset, nt);
      gRNAmeta g;
      g.start_pos = genome_offset_counter - pam.size() - nt;
      g.end_pos = g.start_pos + nt;
      
      auto it = gRNAs.insert(std::make_pair(gRNA, 1));
      if (it.second)
      {
              it.first->second++;
      }
      // std::cout << "Found guide" << std::endl;
      // store metadata
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
  // return 0;
  std::string chr;
  while (1)
  {
    int state = readFastaLine(fp, buf, buf_size - offset, &buf_read);
    // std::cout << "State " << state << std::endl;

    if (state == CHROMOSOME_STATE)
    {
      // read the fasta
      char chromname[1024];
      fgets(chromname, 1024, fp);

      std::size_t pos = chromname.find_first_of(" ");
      chr = chromname;
      if(pos != std::string::npos)
      {
        chr.substr(0,pos);
      }
      
      std::cout << " Processing chromosome: " << chromname << " as " << chr << std::endl;
      
      current_chromosome = buffer;
      buf = buffer;
      pos_counter = 0;
      continue;

    }
    extractgRNA(buffer, buf_read, nt, pam, pos_counter, chr, pos_counter);
    

    // Finally:
    memcpy(buffer, buf + buf_read - offset, offset);
    buf = buffer + offset;
    if (state == 0)
    {
      break;
    }
  }
  
  return 1;
}



int main(int argc, char* argv[])
{
  std::string file = argv[1];
  std::cout << "Reading " << file << std::endl;
  FILE *fp = fopen(file.c_str(), "r");
  if (fp)
  {
    worker(fp,19,"NGG");
  }
  size_t t = 0;
  std::cout << gRNAs.size() << "found" << std::endl;
  for ( auto& gRNA : gRNAs)
  {
    t += gRNA.second;
  }
  std::cout << t << " hits " << std::endl;
  return 0;
}


