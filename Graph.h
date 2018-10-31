#pragma once 
#ifndef GRAPH_H
#define GRAPH_H
#include<string>
#include "Constants.h"

using namespace std;
class Graph {
private:
	double **edges;
	const int V;
public:
	Graph(int rV, string fileName, int seed);
	Graph(Graph& g);
	Graph(int rV);
	int numberOfNodes() { return this->V;}
	void readGraphFromFile(string fileName);
	void writeGraphFile(string fileName);
	~Graph() {
		for (int i = 0; i < V;++i)
			delete[] edges[i];
		delete[] edges;
	}
	double getEdge(int i, int j) {
		if (i >= V || j >= V)
			return INF;
		return edges[i][j];
	}
};
#endif