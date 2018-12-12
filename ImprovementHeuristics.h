#pragma once
#ifndef IH_H
#define IH_H

#include "Graph.h"
#include "Tour.h"
#include "stdafx.h"
#include <ctime>
#include<iostream>
using namespace std;

class ImprovementHeuristics {
public:
	static void TwoOpt(Tour& t, Graph& g);
	static void  ThreeOpt(Tour& t, Graph& g);
};


#endif 
