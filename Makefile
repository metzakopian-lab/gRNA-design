
all:
	g++ -c Trie.cpp
	g++ -c Serialization.cpp
	g++ -O3 -o gRNA -lboost_serialization Serialization.o find-gRNAs.cpp
