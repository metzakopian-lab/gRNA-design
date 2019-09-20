
all:
	g++ -c Trie.cpp
	g++ -c Serialization.cpp
	g++ -O3 -o gRNA -lboost_serialization -fopenmp Serialization.o find-gRNAs.cpp
	g++ -O3 -o build_trie -lboost_serialization Serialization.o Trie.o build-trie.cpp
	g++ -O3 -o gRNA_max -lboost_serialization Serialization.o gRNAs-max.cpp

