


all:
	g++ -c Trie.cpp
	g++ -O3 -o gRNA find-gRNAs.cpp
