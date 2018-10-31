#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H
#include"Graph.h"
#include<vector>
void generateShortestPathInstance(const int NumberOfNode, const int NumberOfInstances);
void readShortestPathInstance(const int NumberOfNode, const int instanceIndex, Graph& g);
#endif