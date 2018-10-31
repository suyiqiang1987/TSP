#include "stdafx.h"
#include "ACS.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include<iomanip>
#include "ConstructiveHeuristics.h"
//#include "ImprovementHeuristics.h"
using namespace std;

ACS::ACS(Graph&g, int rm, double ralpha, double rbeta, double rpho, double rq0 ) : bestTour(g.numberOfNodes()){
	m = rm;
	V = g.numberOfNodes();
	alpha = ralpha;
	beta = rbeta;
	pho = rpho;
	q0 = rq0;
	//Initialize best Tour, tau0, tau and eta
	initialization(g);
}

void ACS::initialization(Graph& g) {
	cout << "...Intilize ACS with Nearest Neighbout Heuristics" << endl;
	Tour nn_t = ConstructiveHeuristics::NNGreedyHeuristic(0, g);
	bestTour.Copy(nn_t);
	tau0 = 1.0 / bestTour.gettotalDistance();
	cout << "...Initial Best Tour:" << bestTour.printTour();
	cout << "...Initial Tau0:" << tau0 << endl;
	cout << "...Initial (i,j):(Tau, Eta):" << endl;
	tau = new double*[V];
	eta = new double*[V];
	for (int i = 0; i < V;i++) {
		cout << "      ";
		tau[i] = new double[V];
		eta[i] = new double[V];
		for (int j = 0; j < V; j++) {
			tau[i][j] = tau0;
			eta[i][j] = 1.0 / g.getEdge(i, j);
			cout << "(" << i << "," << j << "):(" <<setprecision(2) << tau[i][j] << "," <<setprecision(2)<<eta[i][j] << ") ";
		}
		cout << endl;
	}
}