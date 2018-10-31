#include "stdafx.h"
#include "Utilities.h"
#include <fstream>
#include<sstream>
#include <ctime>
#include<string>
#include<iostream>
#include<iomanip>
#include "Graph.h"

using namespace std;

string latlongfile = "LatLongFile.csv";

void generateShortestPathInstance(const int NumberOfNode, const int NumberOfInstances) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Generate " << NumberOfInstances << " Shortest Path Instances with " << NumberOfNode << " Nodes ." << endl;
	for (int i = 0; i < NumberOfInstances;i++) {
		int SeedNumber = i;
		stringstream fileNameStream;
		fileNameStream << (i) << "_" << NumberOfNode  << ".csv";
		string fileName = fileNameStream.str();

		begin = clock();
		cout << "...Generate Instance: Seed = " << i << "; Number of Nodes = " << NumberOfNode << ";" << endl;
		Graph g(NumberOfNode,latlongfile, SeedNumber);
		end = clock();
		timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
		cout << "...Finish Generating Data in " << timeSec << "s for " << i << endl;

		begin = clock();
		g.writeGraphFile(fileName);
		end = clock();
		timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
		cout << "...Finish Writing Data in " << timeSec << "s for " << i << endl;
	}
}

void readShortestPathInstance(const int NumberOfNode, const int instanceIndex, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;

	stringstream fileNameStream;
	fileNameStream << (instanceIndex) << "_" << NumberOfNode <<".csv";
	string fileName = fileNameStream.str();

	begin = clock();
	cout << "Read Instance: Seed = " << instanceIndex << "; Number of Nodes = " << NumberOfNode << endl;
	g.readGraphFromFile(fileName);
	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish Reading Data in " << timeSec << "s for " << instanceIndex << endl;
}