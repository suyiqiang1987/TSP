#pragma once
#ifndef IH_H
#define IH_H

#include "Graph.h"
#include "Tour.h"
#include "stdafx.h"
#include <ctime>
#include<iostream>
using namespace std;

//2-opt General Idea
//Every iteration, the algorithm go through all possibe combination of (i, j), where 0<i<V-1, i+1<=j<V
//Then connect i-1 with j, revere i to j sequence, and connect i to j+1
//Caculate the new cost
//Each iteration update the current tour with best founding 2-opt tour
//The algorithm terminate after fixed number of iteration or after improvement is too small (<0.01)
void TwoOpt(Tour& t, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start 2-Opt..." << endl;
	begin = clock();

	int V = g.numberOfNodes();
	//If there are only two nodes (e.g., origin plus another node), it is no need to use TwoOpt
	if (V < 3) {
		return;
	}

	int iter = 0;
	double improve = 1;
	while (iter < 20 && improve > 0.01) {
		improve = 0;
		Tour bestTour(V);
		bestTour.Copy(t);
		double initial_dist = t.gettotalDistance();
		for (int i = 1; i < V - 1; i++) {
			for (int j = i + 1; j < V; j++) {
				Tour newTour(t);
				newTour.TwoOptSwap(i, j, g);
				if (bestTour.gettotalDistance() > newTour.gettotalDistance()) {
					improve = t.gettotalDistance() - newTour.gettotalDistance();
					bestTour.Copy(newTour);
					//cout << "Swat " << i << " and " << j << " :" << t.printTour();
				}
			}
		}
		if (improve > 0)
			t.Copy(bestTour);
		iter += 1;
		improve /= initial_dist;
	}

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish 2-Opt in " << timeSec << "s." << endl;
}


void ThreeOpt(Tour& t, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start 3-Opt..." << endl;
	begin = clock();

	int V = g.numberOfNodes();
	//If there are only two nodes (e.g., origin plus another node), it is no need to use TwoOpt
	if (V < 3) {
		return;
	}

	//If there are less than 5 nodes, it can only implement 2-opt
	if (V < 5) {
		TwoOpt(t, g);
		return;
	}

	int iter = 0;
	double improve = 1;
	while (iter < 2 && improve > 0.01) {
		improve = 0;
		Tour bestTour(V);
		bestTour.Copy(t);
		double initial_dist = t.gettotalDistance();
		for (int i = 1; i <= V - 4; i++) {
			for (int j = i + 2; j <= V - 2; j++) {
				for (int k = j + 2; k <= V;k++) {
					Tour newTour(t);
					newTour.ThreeOptSwap(i, j, k, g);
					if (bestTour.gettotalDistance() > newTour.gettotalDistance()) {
						improve = t.gettotalDistance() - newTour.gettotalDistance();
						bestTour.Copy(newTour);
						//cout << "Swat " << i << "," << j << " and " << k << ":"<< bestTour.printTour();
					}
				}
			}
		}
		if (improve > 0)
			t.Copy(bestTour);
		iter += 1;
		improve /= initial_dist;
	}

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish 3-Opt in " << timeSec << "s." << endl;
}

#endif 