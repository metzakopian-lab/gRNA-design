#include "Trie.h"

#include <cstring>



Node& Node::expandChildren(const Base& nucl, int id, bool isLast)
{
  // initialize the node
  if(children[nucl] == NULL)
  {
    this->children[nucl] = new Node();

  }
  
  if(isLast)
  {
    this->children[nucl]->incr_node(id);
    this->children[nucl]->isLeaf = true;
  }
  else
  {
    this->children[nucl]->incr_node();
  }

  return *(this->children[nucl]);
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
  memset(children, 0, sizeof(children));
}

