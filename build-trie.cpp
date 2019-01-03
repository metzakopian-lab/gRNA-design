#include <iostream>
#include <map>
#include <vector>
#include <deque>

#include "Trie.h"
#include "Serialization.h"



Node root;
GuideModel gRNAs;

typedef std::map<std::string, int> GuideID;

std::map<std::string, std::vector< std::pair<int,short> > > mismatches;
GuideID guide_uid;


void build_trie(const GuideModel& m, Node& trie_root)
{


  std::cout << "Building trie! " <<std::endl;
  int gid = 0;
  for (auto const& g : m)
  {
    const auto& seq = g.first;
    const auto& meta_vector = g.second;
    guide_uid[g.first] = gid++;
    
    Node* current = &trie_root;
    for(auto const& c : seq)
    {
      // std::cout << (current == &trie_root) << std::endl;
      current = current->expandChildren(map_base(c));
    }
    current->incr_node(guide_uid[g.first]);

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


void mismatch_model(const Node& trie_root, const GuideID& gRNAs, const std::string& guide_str, unsigned short max_mismatches)
{
  std::deque<const Node*> queue;
  Base* enum_base = new Base[guide_str.size()];
  for (auto i = 0; i < guide_str.size(); i++)
  {
    enum_base[i] = map_base(guide_str[i]);
  }
  
  mismatch_rec(queue, trie_root, enum_base, max_mismatches, guide_str.size());
  delete enum_base;
  if(not queue.empty())
  {
    // std::cerr << "Found " << queue.size() << " gRNAs with at least " << max_mismatches << std::endl;
  }
  
}









int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    std::cerr << "Invalid number of arguments" << std::endl;
    std::cout << "Usage:" << std::endl;
    return 1;
  }

  std::string binary_object_filename = argv[1];
  

  modelDeserialize(gRNAs, binary_object_filename);
  build_trie(gRNAs, root);

  // for ( int i = 0 ; i < 4; i++)
  // {
  //   std::cout << root.children[i] << std::endl;
  // }
  int done = 0;
  for(auto const& g : gRNAs)
  {
    mismatch_model(root, guide_uid, g.first, 3);
    
    // std::cout << "." << std::endl;
    if ((++done % 10000) == 0)
      std::cerr << "Processing "<< done << "\r";
  }
  return 0;
}
