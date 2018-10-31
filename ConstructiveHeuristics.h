#pragma once
#ifndef CH_H
#define CH_H

#include "stdafx.h"
#include "Tour.h"
#include "Graph.h"
#include "Constants.h"
#include<ctime>
#include<iostream>
#include <algorithm>    // std::next_permutation, std::sort
using namespace std;

class ConstructiveHeuristics {
public:
	static Tour NNGreedyHeuristic(int startingIndex, Graph& g);
	static Tour Enumeration(int startingIndex, Graph& g);
	static Tour CheapestInseration(int startingIndex, Graph& g);
};


#endif