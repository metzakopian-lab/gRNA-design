#include <vector>
#include <iostream>


enum Base {A, T, G, C};

inline Base map_base(const char& base)
{
  switch (toupper(base))
  {
    case 'A':
      {
        return A;
      }
    case 'T':
      {
        return T;
      }
    case 'G':
      {
        return G;
      }
    case 'C':
      {
        return C;
      }
    default:
      std::cerr << "Unexpected sequence " << base << std::endl;
      return A;
  }
}


class Node
{
  // Dummy class
public:
  int gRNAs;
  bool isLeaf;
  Node* children[4];
  void incr_node(int id);
  void incr_node();
  // for the leaf nodes
  std::vector<int> gRNA_ids;
  
  Node();
  Node* expandChildren(const Base& nucl);
  
  
};



