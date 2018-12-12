#include "stdafx.h"
#include "Graph.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

Graph::Graph(int rV, string fileName, int seed) :V(rV){
	srand(seed);

	edges = new double*[V];
	for (int i = 0; i < V;i++)
		edges[i] = new double[V];

	ifstream file(fileName);
	string line = "";
	// Iterate through each line and split the content using delimetes
	int counter = 0;
	vector<LatLong> latlongs ;
	while (getline(file, line) || counter < V)
	{
		//cout << line << endl;
		istringstream sline(line);
		string s;

		int temp = (rand() % 100);
		if (temp < 20)
			continue;
        
		double latitude;
		double longitude;
		for (int i = 0;i < 2;i++) {
			getline(sline, s, ',');
			if (i == 0) {
				latitude = stod(s);
			}
			else if (i == 1) {
				longitude = stod(s);
				LatLong l(latitude, longitude);
				latlongs.push_back(l);
			}
		}
		counter++;
	}

	for (int i = 0; i < V; i++) {
		//cout << latlongs[i].Latitude << ";" << latlongs[i].Longitude << endl;
		LatLong l1 = latlongs[i];
		for (int j = i; j < V;j++) {
			LatLong l2 = latlongs[j];
			double dist = Constants::haversine(l1, l2);
			this->edges[i][j] = dist;
			if (i != j)
				this->edges[j][i] = dist;
		}
	}


	file.close();
}

Graph::Graph(const Graph&g):V(g.V){
	edges = new double*[V];
	for (int i = 0; i < V;i++) {
		edges[i] = new double[V];
		for (int j = 0; j < V;j++)
			edges[i][j] = g.edges[i][j];
	}
}
Graph::Graph(int rV):V(rV) {
	edges = new double*[V];
	for (int i = 0; i < V;i++)
		edges[i] = new double[V];
}
void Graph::readGraphFromFile(string fileName){
	ifstream file(fileName);
	string line = "";
	// Iterate through each line and split the content using delimeter
	while (getline(file, line))
	{
		//cout << line << endl;
		istringstream sline(line);
		string s;
		int fromIndex;
		int toIndex;
		for (int i = 0;i < 3;i++) {
			getline(sline, s, ',');
			if (i == 0) {
				fromIndex = stoi(s);
			}
			else if (i == 1) {
				toIndex = stoi(s);
			}
			else {
				double d = stod(s);
				//cout << "here" << endl;
				this->edges[fromIndex][toIndex] = d;
			}
		}
	}
	file.close();
}


void Graph::writeGraphFile(string fileName) {
	cout << "...write File " << fileName << endl;
	ofstream myfile;
	myfile.open(fileName);
	for (int i = 0; i < this->V;++i) {
		for (int j = 0; j < this->V;++j) {
			string row_val = to_string(i) + "," + to_string(j) + "," + to_string(this->edges[i][j]) + "\n";
			myfile << row_val;
		}
	}
	myfile.close();
}
