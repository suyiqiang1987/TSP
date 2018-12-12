#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H
#include"Graph.h"
#include<vector>
void generateTSPInstance(const int NumberOfNode, const int NumberOfInstances);
void readTSPInstance(const int NumberOfNode, const int instanceIndex, Graph& g);
#endif
