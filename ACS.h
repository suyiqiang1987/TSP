#pragma once
#ifndef ACS_H
#define ACS_H
#include "Graph.h"
#include "Tour.h"
#include <random>
using namespace std;
typedef pair<double, double> range;
//Ant Colony System is a very effective meta-heuristics for travelling sales man problem
//Ant Colony System is an iteractive algorithm.  At each iteration, a number of artificial ants are considered. Each of them
//builds a solution by walking from node to node on the graph with the constraint of not visiting any vertext that she has already visited in her walk.
//At each stp of the solution construction, an ant selects the following node to be visited according to a stochastic mechanism that is biased by pheromone
class ACS {
private:
	unsigned seed;//ACS ant choose next node to visit according to a random distribution, the seed set the ramdom variable 

	int originIndex; //the node index of starting node
	int numberOfIterations; //Total number of iterations for the ACS
	int m; //Number of artificial ants are considered at each iteration, default = 10
	int V; //Number of nodes
	double alpha;//global pheromone evaporation parameters,default = 0.1
	double beta;//parameter determines the relative importance of pheronmone versus distance, default = 2
	double pho;//local pheromone evaporation parameter,defualt = 0.1
	double q0;//parameter for psudo-random proportional rule, default = 0.9
	double tau0;//global pheromone update parameters, = 1/(V.Lnn), where L is the length of nearest neighbour heuristic results. 
	double **tau;//2-d array store the pheromone of each edge
	double **eta;//2-d array store the inverse of edge length
	int hasNotUpdateBestTourIteration;//Coutner track how many iteration best tour has not been updated
	Tour bestTour; //Best tour found
	void initialization(Graph& g); //Use graph to create eta, and use nearest neighbout heuristic to construct initial solution and update, tau and tau0
	void ConstructAntPath(Graph& g); //Each ant construct a tsp path according the ACS rule
	void localPheromoneUpdate(vector<Tour>& iterTours);//For each ant path, update pheromone to avoid solution converge to a common path to quickly; This is a batch local update (update after each iteration)
	void localPheromoneUpdate(int i, int j);//This is a individual update, update immediately after the edge has been taken;
	void GlobalPheromoneUpdate();//For each iteration, update pheromone to enfornce the good solution
	int chooseNextNode(int currendNodeIndex, const bool Jlist[]); //J is a list of bool value indicate that whether the node has been visited
	double generatePorbabilityDistributionNumerator(int currendNodeIndex, const bool Jlist[], double num[]);//num[] is a list containing of tau * eta^beta
	void generateProbabilityDistribution(double num[], double total, double prob[]); // prob is the cumulative distribution to visit each nodes, the probability to visit node i is equal to prob[i] -prob[i-1]
	double calculateAverageTau(); //Function used to debug the algorithm, the average tau should monotonically decreases 
public:
	ACS(Graph&g, int roriginIndex, int rnumberOfIterations, int rseed, int rm = 10, double ralpha = 0.1, double rbeta = 2, double rpho = 0.1, double rq0 = 0.9);
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
