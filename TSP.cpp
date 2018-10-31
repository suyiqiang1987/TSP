// TSP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Constants.h";
#include "Graph.h"
#include "Tour.h"
#include "Utilities.h"
#include "ConstructiveHeuristics.h"
#include "ImprovementHeuristics.h"
#include "ACS.h"
#include<iostream>
using namespace std;

enum Mode { InstanceGeneration,opt_2,AntColonyOptimization };


int main()
{
	//LatLong l1 = { 36.12, -86.67 };
	//LatLong l2 = { 33.94, -118.4 };
	//cout << "Distance is " << haversine(l1, l2) << endl;
	Graph(50, "LatLongFile.csv",1);
	Mode  m = opt_2;
	int NumberOfNode = 15;
	int NumberOfInstance = 1;
	switch (m)
	{
	case InstanceGeneration:
		generateShortestPathInstance(NumberOfNode, NumberOfInstance);
		break;
	case opt_2:
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readShortestPathInstance(NumberOfNode, i, g);

			Tour nn_t = ConstructiveHeuristics::NNGreedyHeuristic(0, g);
			cout <<"...NN Tour"<<nn_t.printTour();

			//Tour ci_t = ConstructiveHeuristics::CheapestInseration(0, g);
			//cout << "...Cheapest Insertion Tour" << ci_t.printTour();

			Tour two_t = Tour(nn_t);
			TwoOpt(two_t, g);
			cout << "...Two Opt Tour" << two_t.printTour();

			Tour three_t = Tour(nn_t);
			ThreeOpt(three_t, g);
			cout << "...Three Opt Tour" << three_t.printTour();

			//Tour en_t = Enumeration(0, g);
			//cout << "...Enumeration Tour" << en_t.printTour();

			cout << endl;
		}
		break;
	case AntColonyOptimization:
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readShortestPathInstance(NumberOfNode, i, g);

			ACS acs(g);
		}

		break;
	default:
		break;
	}

	char a;
	cin >> a;
    return 0;
}

