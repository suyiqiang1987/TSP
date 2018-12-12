// TSP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Constants.h"
#include "Graph.h"
#include "Tour.h"
#include "Utilities.h"
#include "ConstructiveHeuristics.h"
#include "ImprovementHeuristics.h"
#include "ACS.h"
#include "GENIUS.h"
#include "GENIUS_GPU.h"
#include<iostream>
using namespace std;

enum Mode { InstanceGeneration =0,opt_2=1,AntColonyOptimization=2,GENIUSAlgorithm=3, GENIUSGPUAlgorithm = 4};
//0_48 33551

int main(int argc, char *argv[])
{
	//LatLong l1 = { 36.12, -86.67 };
	//LatLong l2 = { 33.94, -118.4 };
	//cout << "Distance is " << haversine(l1, l2) << endl;
	cout << argv[1] << " " << argv[2] << endl;
	Graph(50, "LatLongFile.csv",1);
	Mode  m = (Mode) stoi(argv[1]);
	int NumberOfNode = stoi(argv[2]);
	int NumberOfInstance = 1;
	int numberOfIterations = 26;
	switch (m)
	{
	case InstanceGeneration:
		generateTSPInstance(NumberOfNode, NumberOfInstance);
		break;
	case opt_2:
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readTSPInstance(NumberOfNode, i, g);

			Tour nn_t = ConstructiveHeuristics::NNGreedyHeuristic(0, g);
			cout <<"...NN Tour"<<nn_t.printTour();

			//Tour ci_t = ConstructiveHeuristics::CheapestInseration(0, g);
			//cout << "...Cheapest Insertion Tour" << ci_t.printTour();

			Tour two_t = Tour(nn_t);
			ImprovementHeuristics::TwoOpt(two_t, g);
			cout << "...Two Opt Tour" << two_t.printTour();

			Tour three_t = Tour(nn_t);
			ImprovementHeuristics::ThreeOpt(three_t, g);
			cout << "...Three Opt Tour" << three_t.printTour();

			//Tour en_t = Enumeration(0, g);
			//cout << "...Enumeration Tour" << en_t.printTour();

			cout << endl;
		}
		break;
	case AntColonyOptimization:
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readTSPInstance(NumberOfNode, i, g);

			//Tour en_t = ConstructiveHeuristics::Enumeration(0, g);
			//cout << "...Enumeration Tour" << en_t.printTour();

			ACS acs(g,0, numberOfIterations,1);

			cout << endl<<endl;
		}

		break;

	case GENIUSAlgorithm:
		cout << "Opt 3: Use GENI-US Algorithm" << endl;
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readTSPInstance(NumberOfNode, i, g);

			//Tour en_t = ConstructiveHeuristics::Enumeration(0, g);
			//cout << "...Enumeration Tour" << en_t.printTour();

			GENIUS genius(g, 0);
			

		}
		break;

	case GENIUSGPUAlgorithm:
		cout << "Opt 3: Use GENI-US Algorithm" << endl;
		for (int i = 0; i < NumberOfInstance;i++) {
			Graph g(NumberOfNode);
			readTSPInstance(NumberOfNode, i, g);

			//Tour en_t = ConstructiveHeuristics::Enumeration(0, g);
			//cout << "...Enumeration Tour" << en_t.printTour();

			GENIUS_GPU genius(g, 0);
		}
		break;
	default:
		break;
	}

	char a;
	cin >> a;
    return 0;
}

