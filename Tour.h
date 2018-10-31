#pragma once
#ifndef TOUR_H
#define TOUR_H
#include "Graph.h"
#include <string>
class Tour {
private:
	int V; //number of nodes
	int* seq;//Sequence node of the tour, the first node is origin, the default destination will be the same as origin
	bool* seenList;//List of boolean value to indicate whether the node has been added to the tour
	double* weight;//weight[0] is the distance between seq[0] and seq[1],..., weight[V-1] is the distance between seq[V-1] and destination
	int lastIndex; //This value is used to store information of constructive heuristic, which the number of nodes that have been added to the tour
	double totalDistance; //The total distance of tour of nodes already added to the tour
public:
	Tour(int rV);
	Tour(Tour& t);
	void Copy(Tour& t);
	int getSeq(int index) { return seq[index];}
	bool getHasBeenSeen(int index) {return seenList[index];}
	double gettotalDistance() { return totalDistance; }
	int getLastIndex() { return lastIndex; }
	double calculateAddNodeAt(int i, int location, Graph& g);//Chcek the cost of adding node i at location, and return the total cost;
	bool addNode(int i, int location, Graph& g);//Add the node i into the location, and update all the parameters; If the node cannot be added to the tour, then return False;
	bool canAdd(int i, int location);
	void TwoOptSwap(int i, int j, Graph& g);
	void ThreeOptSwap(int i, int j, int k,Graph& g);
	string printTour();
	~Tour() {
		delete[] seq;
		delete[] weight;
		delete[] seenList;
	}
};
#endif