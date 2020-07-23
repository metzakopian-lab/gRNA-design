
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <ctype.h>
#include <cstdio>
#include <cstring>
#include <omp.h>
#include "Model.h"
#include "Serialization.h"

#define BUFFER_SIZE 120000

#define N_STATE 1
#define SEQ_STATE 1 
#define CHROMOSOME_STATE 2

GuideModel gRNAs;
char inv_nucl[256];



int readFastaLine(std::ifstream& fasta, std::string& buffer)
{
  
  unsigned int chunk = 100 * (2 << 20);
  buffer.reserve(chunk);
  unsigned int bytes_left = chunk;
  
  size_t i;
  
  if (fasta.peek() == '>')
  {
    
    std::getline(fasta, buffer);
    return CHROMOSOME_STATE;
  }
 
  while (!fasta.eof() && (fasta.peek() != '>'))
  {
    std::string tmp;
    tmp.reserve(256);
    std::getline(fasta, tmp);
    // std::cout << tmp << std::endl;
    // Remove character?
    // std::cerr << "Removed " << tmp[tmp.size() -1] << " ";
    // tmp.pop_back();
    
    buffer += tmp;
    bytes_left -= tmp.size();
    if (bytes_left  <  tmp.size() * 100 )
    {
      buffer.reserve(buffer.size() + chunk);
      bytes_left = chunk;
    }
  }
  return SEQ_STATE;
}


// Extracting chromosomes from 
void extractgRNA(const std::string& buf,  int nt, const std::string& pam, const std::string& chr_name, unsigned int& guide_id)
{
  const unsigned int total_seq_length = nt + pam.size();

  std::string inv_pam = pam;
  for (int i = 0; i < pam.size(); ++i)
  {
    char cc = pam[i];
    inv_pam[inv_pam.size() - 1 - i] = inv_nucl[cc];
  }
  
  for (auto i = total_seq_length; i < buf.size(); ++i)
  {

    bool inv_state = true;
    bool for_state = true;
    
    if (buf[i] == 'N')
    {
      
      i += total_seq_length - 1;
      
      continue;
    }

    // Check sense and antisense strand
    for ( unsigned int p = 0; p < pam.size(); ++p)
    {
      
      if (pam[p] !=  'N' && (pam[p] != buf[i - pam.size() + p]))
      {
        for_state = false;
      }

      if (inv_pam[p] != 'N' && (inv_pam[p] != buf[i - total_seq_length + p]))
      {
        inv_state = false;
      }
    }

    if (inv_state || for_state) 
    {
      GuideMeta g;
      bool invalid = false;
      auto gRNA = buf.substr(i - total_seq_length, total_seq_length);

      // Filter N's (optimized)
      for (auto const& c : gRNA)
      {
        if (c == 'N')
        {
          invalid = true;
          break;
        }

      }
      if (invalid)
      {
        continue;
      }
      
      for ( unsigned int p = 0; p < pam.size(); ++p)
      {
        if(pam[p] == 'N' && for_state)
        {
           gRNA[gRNA.size() - pam.size() + p] = 'N' ;
          
        }
        if (inv_pam[p] == 'N' && inv_state)
        {
          gRNA[p] = 'N';

        }
      }
      g.chromosome = chr_name;
      g.end_pos = i;
      g.start_pos = g.end_pos - total_seq_length;
      // std::cout << ">" << g.id <<"," << (inv_state ? "inv" : "for")  << g.chromosome << ":" << g.start_pos << "-" << g.end_pos << std::endl;
      // std::cout << gRNA << std::endl;
      
      
      #pragma omp critical
      {
        g.id = guide_id++;
        auto it = gRNAs.insert(std::make_pair(gRNA, std::vector<GuideMeta>()));
        if (it.second)
        {
          it.first->second.push_back(g);
        }
      }
      

      
    }
  }
}

int worker(std::ifstream& fasta, const unsigned char nt, const std::string pam)
{
  
  unsigned int gRNA_id = 0;
  
  std::vector <std::pair<std::string, std::string>> chromosomes;
  std::cerr << "Reading Fasta";
  while (!fasta.eof())
  {
    std::pair<std::string, std::string> chr;
    int state1 = readFastaLine(fasta, chr.first);
    int state2 = readFastaLine(fasta, chr.second);
    
    if (not (state1 == CHROMOSOME_STATE && state2 == SEQ_STATE))
    {
      std::cerr << "Something went wrong" << std::endl;
      return 1;
    }
    std::string chr_name;
    std::size_t pos = chr.first.find_first_of(" \t\n");
    if (pos != std::string::npos)
    {
      chr_name = chr.first.substr(1, pos - 1);
    }
    else
    {
      chr_name = chr.first.substr(1);
    }
    chr.first = chr_name;
    // std::cerr << "Read chromosome: " << chr.first << " as " << chr_name << std::endl; 
    chromosomes.push_back(chr);
  }
  
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < chromosomes.size(); i++)
  {
    std::cerr << "Processing chromosome: " << chromosomes[i].first << " " << chromosomes[i].second.size() / 1000000 << "mBases" << std::endl; 
    extractgRNA(chromosomes[i].second, nt, pam, chromosomes[i].first, gRNA_id);
    std::cerr << "Finished chromosome: " << chromosomes[i].first << std::endl; 
    chromosomes[i].second.clear();
  }

  
  

  #pragma omp
  return 0;
}



int main(int argc, char* argv[])
{

  inv_nucl['A'] = 'T';
  inv_nucl['T'] = 'A';
  inv_nucl['G'] = 'C';
  inv_nucl['C'] = 'G';
  inv_nucl['N'] = 'N';


  if(argc != 3)
  {
    std::cerr << "Invalid number of arguments" << std::endl;
    std::cerr << "Usage: " << argv[0] << "<fasta file> <fasta-output-file>" << std::endl;
    return 1;
  }

  std::string file = argv[1];
//  std::string gRNA_file = argv[2];
  std::string gRNAs_fasta = argv[2];
  
  std::cerr << "Input file " << file << std::endl;
  std::ifstream  fp(file.c_str());
  
  if (fp)
  {
    worker(fp, 19, "NGG");
  }
  else
  {
    std::cerr << file << " not found" << std::endl;
    return 1;
  }
  
  //modelSerialize(gRNAs, gRNA_file);


  std::ofstream out(gRNAs_fasta.c_str());
  if (out)
  {
    for(auto const & _g : gRNAs)
    {
      auto g = _g.second[0];
      auto gRNA = _g.first;
      std::cout << ">" << g.id <<"," << g.chromosome << ":" << g.start_pos << "-" << g.end_pos << std::endl;
      std::cout << gRNA << std::endl;
    }

  }
  else
  {
    std::cerr << "Could not open to" << gRNAs_fasta.c_str() <<std::endl;
  }
 
  return 0;
}


