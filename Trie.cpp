#include <cstring>

#include "Trie.h"





Node* Node::expandChildren(const Base& nucl)
{
  // initialize the node
  if(children[nucl] == NULL)
  {
    this->children[nucl] = new Node();

  }

  this->children[nucl]->incr_node();
  
  return this->children[nucl];
}

void Node::incr_node()
{ 
  this->gRNAs++; 
}

void Node::incr_node(int id)
{
  gRNA_ids.push_back(id);
}


Node::Node() : isLeaf(false)
{
 
  for (auto i = 0; i < 4; i++){
    children[i] = NULL;
  }
}

