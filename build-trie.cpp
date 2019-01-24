#include <iostream>
#include <map>
#include <vector>
#include <deque>
#include <fstream>

#include "Trie.h"
#include "Serialization.h"



Node root;
GuideModel gRNAs;

typedef std::map<std::string, int> GuideID;
typedef std::map<int,std::string> InvGuideID;

std::map<int, std::vector< std::pair<int, int> > > mismatches;
GuideID guide_uid;
InvGuideID inv_guide_uid;


void build_trie(const GuideModel& m, Node& trie_root)
{


  std::cout << "Building trie! " <<std::endl;
  int gid = 0;
  for (auto const& g : m)
  {
    const auto& seq = g.first;
    int guide_id = gid++;
    guide_uid[g.first] = guide_id;
    inv_guide_uid[guide_id] = g.first;
    
    Node* current = &trie_root;
    for(auto const& c : seq)
    {
      current = current->expandChildren(map_base(c)); 
    }
    current->incr_node(guide_id);

    current->isLeaf = true;
    
  }
}


void mismatch_rec(std::deque<const Node*>& queue,
                                     const Node& trie_root,
                                     const Base* gRNA,
                                     const short mismatches_left,
                                     unsigned short nt_left)
{
  

  if (trie_root.isLeaf && 
      nt_left == 0 &&
      mismatches_left >= 0)
  {
    queue.push_back(&trie_root);
  }

  for (short i = 0; i < 4; ++i)
  {

    if(trie_root.children[i] == NULL)
    {
      continue;
    }
    // std::cout << nt_left << std::endl;
    const Node& v = *(trie_root.children[i]);
    if(i == gRNA[0])
    {
      mismatch_rec(queue, v, gRNA + 1, mismatches_left, nt_left - 1);
    }
    else
    {
      if(mismatches_left == 0)
      {
        // criteria unmet
        continue;
      }
      else
      {
        mismatch_rec(queue, v, gRNA + 1, mismatches_left - 1, nt_left - 1);
      } 
    }
  }
}


void mismatch_model(const Node& trie_root, int gRNAid, const std::string& guide_str, unsigned short max_mismatches)
{
  std::deque<const Node*> queue;
  Base* enum_base = new Base[guide_str.size()];
  for (auto i = 0; i < guide_str.size(); i++)
  {
    enum_base[i] = map_base(guide_str[i]);
  }
  
  mismatch_rec(queue, trie_root, enum_base, max_mismatches, guide_str.size());
  delete enum_base;
  
  mismatches[gRNAid] = std::vector<std::pair<int,int>>();
  auto& g_mis = mismatches[gRNAid];
  for(auto const& g : queue)
  {
    int score = 0;
    const std::string& mismatch = inv_guide_uid[g->gRNA_ids[0]];
    for (auto i = 0; i < guide_str.size(); i++)
    {
      score |= ((guide_str[i] == mismatch[i]) << i);
    }
    g_mis.push_back(std::make_pair(g->gRNA_ids[0], score));
  } 
}




void export_gIDs(const std::string filename)
{
  std::ofstream out(filename.c_str());
  out << "gRNAid,sequence"<< std::endl;
  for (auto const& g : guide_uid)
  {
    out << g.first << "," << g.second << std::endl;
  }
  out.close();
}

void export_mismatches(const std::string filename)
{
  std::ofstream out(filename.c_str());
  out << "gID1,gID2,score" << std::endl;
  for (auto const& g : mismatches)
  {
    for( auto const& m : g.second)
    {
      out << g.first << "," << m.first << "," << m.second << std::endl;
    }
  }
  out.close();
}





int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    std::cerr << "Invalid number of arguments"<< std::endl;
    std::cout << "Usage: " <<argv[0] <<" <model-file> <gRNA-ids-csv> <mismatches-csv> " << std::endl;
    return 1;
  }

  std::string binary_object_filename = argv[1];
  
  modelDeserialize(gRNAs, binary_object_filename);
  build_trie(gRNAs, root);
  int done = 0;
  for(auto const& g : gRNAs)
  {
    mismatch_model(root, guide_uid[g.first], g.first, 3);
    if ((++done % 10000) == 0)
      std::cerr << "Processing "<< done << "\r";
  }
  std::cerr << "Exporting data" << std::endl;
  
  export_gIDs(argv[2]);
  export_mismatches(argv[3]);
  
  
  return 0;
}
