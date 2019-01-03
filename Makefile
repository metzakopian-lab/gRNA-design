


all:
	g++ -c Trie.cpp
	g++ -O3 -o gRNA -lboost_serialization find-gRNAs.cpp
