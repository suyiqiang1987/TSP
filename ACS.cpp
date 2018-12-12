#include "stdafx.h"
#include "ACS.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include<list>
#include<iomanip>
#include "ConstructiveHeuristics.h"
#include "ImprovementHeuristics.h"
//#include "ImprovementHeuristics.h"

ACS::ACS(Graph&g, int roriginIndex,  int rnumberOfIterations, int rseed, int rm, double ralpha, double rbeta, double rpho, double rq0 ) : bestTour(g.numberOfNodes()){
	clock_t begin;
	clock_t end;
	double timeSec;
	begin = clock();
	
	hasNotUpdateBestTourIteration = 0;
	originIndex = roriginIndex;
	seed = rseed;
	srand(seed);

	numberOfIterations = rnumberOfIterations;
	m = rm;
	V = g.numberOfNodes();
	alpha = ralpha;
	beta = rbeta;
	pho = rpho;
	q0 = rq0;
	cout << "Alpha = " << alpha << "; Beta = " << beta << "; pho = " << pho << "; q0 =" << q0 << endl;
	//Initialize best Tour, tau0, tau and eta
	initialization(g);
	double lastTauAverage = calculateAverageTau();
	double lastIterationBestTourDist = bestTour.gettotalDistance();
	for (int i = 0; i < numberOfIterations; i++) {
		//cout << "Iteration " << i << ":";
		ConstructAntPath(g);
		ImprovementHeuristics::ThreeOpt(bestTour,g);
		GlobalPheromoneUpdate();

		double curTauAverage = calculateAverageTau();
		//double diff = 1 - curTauAverage / lastTauAverage;
		//cout << "..Iteration " << i << "'s avg tau = " << curTauAverage << " which is different from last iteration by " <<setprecision(2) << diff <<" and best tour found so far :" << bestTour.printTour();
		//cout <<  "..Iteration " << i << "'s avg tau = " << curTauAverage << "; Best Tour Lenght = " << bestTour.gettotalDistance() << endl;
		lastTauAverage = curTauAverage;
		
		if (lastIterationBestTourDist - bestTour.gettotalDistance()>0)
			hasNotUpdateBestTourIteration = 0;
		else
			hasNotUpdateBestTourIteration++;
		lastIterationBestTourDist = bestTour.gettotalDistance();

		//if (hasNotUpdateBestTourIteration>10)
			//break;
	}
	cout << "ACS Tour:" << bestTour.printTour();
	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish ACS in " << timeSec << endl;
}

void ACS::initialization(Graph& g) {
	cout << "...Intilize ACS with Nearest Neighbout Heuristics" << endl;
	Tour nn_t = ConstructiveHeuristics::NNGreedyHeuristic(0, g);
	bestTour.Copy(nn_t);
	tau0 = 1.0 / (V*bestTour.gettotalDistance()); // tau_0 = 1/(Number of Nodes * NN Tour Length), this is according to http://people.idsia.ch/~luca/acs-ec97.pdf
	cout << "...Initial Best Tour:" << bestTour.printTour();
	cout << "...Initial Tau0:" << tau0 << endl;
	//cout << "...Initial (i,j):(Tau, Eta):" << endl;
	tau = new double*[V];
	eta = new double*[V];
	for (int i = 0; i < V;i++) {
		//cout << "      ";
		tau[i] = new double[V];
		eta[i] = new double[V];
		for (int j = 0; j < V; j++) {
			tau[i][j] = tau0;
			if (g.getEdge(i, j)>0)
				eta[i][j] = 1.0 / g.getEdge(i, j);
			else
				eta[i][j] = Constants::INF;
			//cout << "(" << i << "," << j << "):(" <<setprecision(2) << tau[i][j] << "," <<setprecision(2)<<eta[i][j] << ") ";
		}
		//cout << endl;
	}
}

void ACS::GlobalPheromoneUpdate() {
	//Every edge's tau evaporate by pho
	for (int i = 0; i < V;i++) {
		for (int j = 0; j < V;j++)
			this->tau[i][j] = this->tau[i][j] * (1 - alpha);
	}
	//The edge along best tour is increased by alpha * (1/best tour length)
	for (int i = 0; i < V; i++) {
		int this_node = bestTour.getSeq(i);
		int next_node = bestTour.getSeq((i + 1) % V);
		this->tau[this_node][next_node] = alpha * (1 / bestTour.gettotalDistance()) ;
	}
}

void ACS::localPheromoneUpdate(vector<Tour>& iterTours) {//For each ant path, update pheromone to avoid solution converge to a common path to quickly; This is a batch local update (update after each iteration)
	vector<Tour>::iterator i;
	for (i = iterTours.begin(); i != iterTours.end();i++) {
		for (int j = 0; j < V;j++) {
			int k = (j + 1) % V;
			int this_node = (*i).getSeq(j);
			int next_node = (*i).getSeq(k);
			localPheromoneUpdate(this_node, next_node);
		}
	}
}


void ACS::ConstructAntPath(Graph& g) {
	vector<Tour> iterTours;
	for (int i = 0; i < m; i++) {
		int currendNodeIndex = rand() % (V - 1) + 1;
		while(currendNodeIndex==originIndex)
			currendNodeIndex = rand() % (V - 1) + 1;
		//cout << "m = " << i << " 1st node :" << currendNodeIndex <<endl;
		bool* Jlist = new bool[this->V];
		for (int k = 0; k < this->V;k++)
			Jlist[k] = false;

		Tour t(this->V);
		t.addNode(originIndex, 0, g);
		Jlist[originIndex] = true;
		t.addNode(currendNodeIndex,1,g);
		Jlist[currendNodeIndex] = true;
		for (int j = 2; j < V; j++) {
			currendNodeIndex = chooseNextNode(currendNodeIndex, Jlist);
			t.addNode(currendNodeIndex, j, g);
			Jlist[currendNodeIndex] = true;
			localPheromoneUpdate(t.getSeq(j - 1), t.getSeq(j));
		}
		if (bestTour.gettotalDistance() > t.gettotalDistance()) {
			bestTour.Copy(t);
			//cout << t.printTour();
		}
		iterTours.push_back(t);

		delete[] Jlist;
	}


}

void ACS::localPheromoneUpdate(int i, int j) {//This is a individual update, update immediately after the edge has been taken;
	this->tau[i][j] = (1 - pho) * this->tau[i][j] + pho * tau0;
}


int ACS::chooseNextNode(int currendNodeIndex, const bool Jlist[]) {
	int s=-1;
	double *num = new double[this->V];
	double total = generatePorbabilityDistributionNumerator(currendNodeIndex, Jlist, num);
	double q = (rand()%10000)/10000.0;//distribution(generator);
	if (q <= q0) {
		s = distance(num, max_element(num, num + V));//this statement return the index of the max element in the list
	}
	else {
		double p = (rand() % 10000) / 10000.0;
		double *prob = new double[this->V];
		generateProbabilityDistribution(num, total, prob);
		//cout << p << endl;
		if (p < prob[0])
			s = 0;
		else {
			for (int i = 1; i < V;i++) {
				//cout << "..." << prob[i];
				//if a node has been visited then prob[i] = prob[i-1], then it will not go to the if condition below
				if (p >= prob[i - 1] && p < prob[i]) {
					s = i;
					break;
				}
			}
		}

		delete[] prob;
	}

	
	delete [] num;
	return s;
}

double ACS::generatePorbabilityDistributionNumerator(int currendNodeIndex, const bool Jlist[], double num[]) {
	double total = 0.0;
	for (int i = 0; i < V;i++) {
		if (Jlist[i] == true) {
			num[i] = 0;
		}
		else {
			num[i] = tau[currendNodeIndex][i] * pow(eta[currendNodeIndex][i], beta);
			//cout <<"...." << num[i] << endl;
			total += num[i];
		}
	}
	//cout << total << endl;
	return total;
}
//Current Node Index is the node's index, not it is sequence index
void ACS::generateProbabilityDistribution(double num[], double total, double prob[]) {
	prob[0] = num[0] / total;
	for (int i = 1; i < V; i++) {
		prob[i] = prob[i - 1] + num[i] / total;
	}
}

double ACS::calculateAverageTau() {
	double avg = 0.0;
	for (int i = 0; i < V;i++) {
		for (int j = 0; j < V;j++) {
			avg += this->tau[i][j];
		}
	}
	avg = avg / (V*V);
	return avg;
}
