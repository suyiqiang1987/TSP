#pragma once
#ifndef GENIUS_H
#define GENIUS_H
#include "Graph.h"
#include "Tour.h"
#include <vector>
#include <random>
using namespace std;
//GENIUS (Generalized Insertion - Unstring - String) algorithm is an effective algorithm to solve TSP Problem. 
class GENIUS {
private:
	int V; //Number of nodes in the graph;
	int p; //Neighbourhood size of searching
	int originIndex; //the node index of starting node
	int totalNumberOfiterations;
	int ** pNeighbourhood;
	double **edges;
	vector<int> unsolvedNode;
	Tour bestTour;
	void updatePNeighbourhood(int newIndex, Graph& g);//When a new node insert to the tour, need to update the pNeighbourhood array;
	bool isInPNeighbourHood(int i, int j);
	void reverse(int i, int j, int* temp_array, int& counter, Tour t);
	void reverse(int i, int j, int* temp_array, int& counter, int* t, int lastIndex);
	void addBetween(int i, int j, int* temp_array, int& counter, Tour t);
	void addBetween(int i, int j, int* temp_array, int& counter, int *t, int lastIndex);
	void Initialization(Graph& g);//According to paper, it needs to initalize the best tour by arbitrarily selecting three vetices. 

	void GENI(Graph& g); //Constructive phase of GENIUS
	void typeIInsertion(Graph& g, int iter, int nodeV, Tour& typeIBestTour, Tour& t);

	void typeIIInsertion(Graph& g, int iter, int nodeV, Tour& typeIIBestTour, Tour& t);

	void US(Graph& g);//Improve Phase of GENIUS
	void typeIUnstring(Graph& g, int i, Tour& typeIUNstringBestTour, Tour& t);
	void typeIIUnstring(Graph& g, int i, Tour& typeIIUNstringBestTour, Tour& t);

public:
	GENIUS(Graph& g, int roriginIndex, int rp = 3, int rtotalNumberOfiterations=200);
	
	~GENIUS() {
		for (int i = 0; i < this->V; i++)
			delete[] pNeighbourhood[i];
		delete[] pNeighbourhood;

		for (int i = 0; i < V;++i)
			delete[] edges[i];
		delete[] edges;
	}
};

#endif
