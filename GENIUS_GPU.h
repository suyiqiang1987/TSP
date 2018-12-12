#pragma once
#ifndef GENIUSGPU_H
#define GENIUSGPU_H
#include "Graph.h"
#include "Tour.h"
#include <vector>
#include <random>
using namespace std;
//GENIUS (Generalized Insertion - Unstring - String) algorithm is an effective algorithm to solve TSP Problem. 
class GENIUS_GPU {
private:
	int V; //Number of nodes in the graph;
	int p; //Neighbourhood size of searching
	int originIndex; //the node index of starting node
	int totalNumberOfiterations;
	int**  pNeighbourhood;
	vector<int> unsolvedNode;
	Tour bestTour;
	void updatePNeighbourhood(int newIndex, Graph& g);//When a new node insert to the tour, need to update the pNeighbourhood array;
	bool isInPNeighbourHood(int i, int j, int neighboursize,int** __restrict pNeighbourhood);
	void reverse(int i, int j, int* temp_array, int& counter, Tour t);
	void reverse(int i, int j, int* temp_array, int& counter, int* t, int lastIndex);
	double reverseDist(int i, int j, int *__restrict t, int lastIndex, double *__restrict*__restrict e);
	void addBetween(int i, int j, int* temp_array, int& counter, Tour t);
	void addBetween(int i, int j, int* temp_array, int& counter, int *t, int lastIndex);
	double addBetweenDist(int i, int j, int *__restrict t, int lastIndex, double *__restrict*__restrict e);
	void Initialization(Graph& g);//According to paper, it needs to initalize the best tour by arbitrarily selecting three vetices. 

	void GENI(Graph& g, double *__restrict*__restrict e, int**  pNeighbourhood); //Constructive phase of GENIUS
	void typeIInsertion(Graph& g, int iter, int nodeV, Tour& typeIBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhoode);
	void typeIInsertionGPU(Graph& g, int iter, int nodeV, Tour& typeIBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhoode);

	void typeIIInsertion(Graph& g, int iter, int nodeV, Tour& typeIIBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhood);

	void US(Graph& g, double *__restrict*__restrict e, int**  pNeighbourhood);//Improve Phase of GENIUS
	void typeIUnstring(Graph& g, int i, Tour& typeIUNstringBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhood);
	void typeIIUnstring(Graph& g, int i, Tour& typeIIUNstringBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhood);

public:
	GENIUS_GPU(Graph& g, int roriginIndex, int rp = 3, int rtotalNumberOfiterations=200);
	
	~GENIUS_GPU() {
		for (int i = 0; i < this->V; i++)
			delete[] pNeighbourhood[i];
		delete[] pNeighbourhood;


	}
};

#endif
