#pragma once
#ifndef ACS_H
#define ACS_H
#include "Graph.h"
#include "Tour.h"
//Ant Colony System is a very effective meta-heuristics for travelling sales man problem
//Ant Colony System is an iteractive algorithm.  At each iteration, a number of artificial ants are considered. Each of them
//builds a solution by walking from node to node on the graph with the constraint of not visiting any vertext that she has already visited in her walk.
//At each stp of the solution construction, an ant selects the following node to be visited according to a stochastic mechanism that is biased by pheromone
class ACS {
private:
	int m; //Number of artificial ants are considered at each iteration, default = 10
	int V; //Number of nodes
	double alpha;//global pheromone evaporation parameters,default = 0.1
	double beta;//parameter determines the relative importance of pheronmone versus distance, default = 2
	double pho;//local pheromone evaporation parameter,defualt = 0.1
	double q0;//parameter for psudo-random proportional rule, default = 0.9
	double tau0;//global pheromone update parameters, = 1/(V.Lnn), where L is the length of nearest neighbour heuristic results. 
	double **tau;//2-d array store the pheromone of each edge
	double **eta;//2-d array store the inverse of edge length
	Tour bestTour; //Best tour found
	void initialization(Graph& g); //Use graph to create eta, and use nearest neighbout heuristic to construct initial solution and update, tau and tau0
	void ConstructAntPath(Tour& t); //Each ant construct a tsp path according the ACS rule
	void localPheromoneUpdate();//For each ant path, update pheromone to avoid solution converge to a common path to quickly
	void GlobalPheromoneUpdate();//For each iteration, update pheromone to enfornce the good solution
public:
	ACS(Graph&g, int rm = 10, double ralpha = 0.1, double rbeta = 2, double rpho = 0.1, double rq0 = 0.9);
	~ACS() {
		for (int i = 0; i < V; i++) {
			delete[] tau[i];
			delete[] eta[i];
		}
		delete[] tau;
		delete[] eta;
	}
};
#endif