#include "stdafx.h"
#include"Tour.h"
#include "Constants.h"
#include<iostream>
#include<string>
#include<math.h>
using namespace std;

Tour::Tour() :V(100) {
	seq = new int[this->V];
	weight = new double[this->V];
	seenList = new bool[this->V];
	for (int i = 0; i < V;++i)
		seenList[i] = false;
	lastIndex = -1;
	totalDistance = 0;
}

Tour::Tour(int rV) :V(rV) {
	seq = new int[this->V];
	weight = new double[this->V];
	seenList = new bool[this->V];
	for (int i = 0; i < V;++i)
		seenList[i] = false;
	lastIndex = -1;
	totalDistance = 0;
}
Tour::Tour(const Tour& t) : V(t.V) {
	seq = new int[this->V];
	weight = new double[this->V];
	seenList = new bool[this->V];
	lastIndex = t.lastIndex;
	totalDistance = t.totalDistance;
	for (int i = 0; i < V; i++) {
		seq[i] = t.seq[i];
		weight[i] = t.weight[i];
		seenList[i] = t.seenList[i];
	}
}

void Tour::Copy(const Tour& t) {
	lastIndex = t.lastIndex;
	totalDistance = t.totalDistance;
	for (int i = 0; i < V; i++) {
		seq[i] = t.seq[i];
		weight[i] = t.weight[i];
		seenList[i] = t.seenList[i];
	}
}

double Tour::calculateAddNodeAt(int i, int location, Graph& g)
{
	bool canAddFlag = canAdd(i, location);
	double val = this->totalDistance;

	if (canAddFlag == false)
		return Constants::INF;
	if (location == 0) {
		weight[0] = 0;
	}
	else if (location == lastIndex + 1) {//If the location is after the last node location, then subtract the distance from
		//last node back to origin, and add distance from last node to node i and distance from node i to origin
		val -= weight[lastIndex];
		val += g.getEdge(seq[lastIndex], i);
		val += g.getEdge(i, seq[0]);//If add i to end of tour, then the weight is the distance back to origin
	}
	else {//Then subtract the distance from location to location + 1, and add distance from location to i, and distance from i to location + 1
		val -= weight[location-1];
		val += g.getEdge(seq[location-1], i);
		val += g.getEdge(i, seq[location]);//If add i to end of tour, then the weight is the distance back to origin
	}

	return val;
}

bool Tour::addNode(int i, int location, Graph& g) {
	//cout<<"Tried to add " << i << " at " << location << endl;
	bool canAddFlag = canAdd(i,location);
	if (canAddFlag == false)
		return false;

	//cout << "Success to add " << i << " at " << location << endl;
	seenList[i] = true;
	for (int index = lastIndex; index >= location;index--) {
		seq[index + 1] = seq[index]; //shift all the value of seq right by one location
		weight[index + 1] = weight[index];//shift all the value of weight right by one location
	}

	seq[location] = i;
	if (location == 0) {
		weight[0] = 0;
	}
	else if (location == lastIndex + 1) {
		weight[location - 1] = g.getEdge(seq[location - 1], i);
		weight[location] = g.getEdge(i, seq[0]);//If add i to end of tour, then the weight is the distance back to origin
	}
	else {
		weight[location - 1] = g.getEdge(seq[location - 1], i);
		weight[location] = g.getEdge(i, seq[location+1]);//If add i to end of tour, then the weight is the distance back to origin
	}
	lastIndex = lastIndex + 1;

	totalDistance = 0;
	for (int index = 0; index <= lastIndex;index++)
		totalDistance += weight[index];

	return canAddFlag;
}

bool Tour::canAdd(int i, int location) {
	if (i >= V || i < 0) {//If i is out of range, or i has been added to the tour, 
		cout << "Node Index" << i << " has to be >=0 and <= " << V << endl;
		return false;
	}

	if (location > lastIndex + 1 || location < 0) {
		cout << "Tour Location " << location << " has to be >= 0 and <=" << (lastIndex + 1) << endl;
		return false;
	}

	if (seenList[i] == true) {
		cout << i << " has already been added to tour" << endl;
		return false;
	}

	if (location == 0 && lastIndex != -1) {
		cout << "Origin " << seq[0] << " has been added, and could not be changed to " << i << endl;
		return false;
	}

	return true;
}

void Tour::TwoOptSwap(int i, int j, Graph& g) {
	//j must be greater than i+1 to implement the 2-opt
	if (i + 1 > j)
		return;

	//Cannt 2-opt origin
	if (i < 1 || i>V-1)
		return;

	//Index out of range;
	if (j >= V)
		return;

	/*int lengthOfReversalList = j - i;
	int* reveralList = new int[lengthOfReversalList+1];
	for (int k = j; k>= i;k--) 
		reveralList[j-k] = seq[k];

	for (int k = i; k <= j;k++) 
		seq[k] = reveralList[k - i];*/

	rvereseArray(seq, i, j);
	
	for (int k = i - 1; k <= j;k++) {
		totalDistance -= weight[k];
		weight[k] = g.getEdge(seq[k], seq[(k+1) % V]);
		totalDistance += weight[k];
	}

	//delete [] reveralList;

}

void Tour::ThreeOptSwap(int i, int j, int k, Graph& g) {
	//If the below constriant is violated, then the index will be out of range
	if (i + 2 > j || j + 2 > k)
		return;
	//Cannot opt origin and destination
	if (i < 1 || i > V-4)
		return;
	

	int temp_arr[4] = { i,j-1,j,k-1};
	Tour bestTour(this->V);
	for (int per_index = 0; per_index < 7; per_index++) {
		Tour temp(this->V);
		int counter = 0;
		for (int seq_index = 0; seq_index <= i - 1;seq_index++) 
			temp.addNode(seq[seq_index], counter++, g);
		
		if (temp_arr[Constants::threeOptPermutation[per_index][0]] < temp_arr[Constants::threeOptPermutation[per_index][1]]) {
			for (int seq_index = temp_arr[Constants::threeOptPermutation[per_index][0]];seq_index <= temp_arr[Constants::threeOptPermutation[per_index][1]];seq_index++)
				temp.addNode(seq[seq_index], counter++, g);
		}
		else {
			for (int seq_index = temp_arr[Constants::threeOptPermutation[per_index][0]];seq_index >= temp_arr[Constants::threeOptPermutation[per_index][1]];seq_index--)
				temp.addNode(seq[seq_index], counter++, g);
		}

		if (temp_arr[Constants::threeOptPermutation[per_index][2]] < temp_arr[Constants::threeOptPermutation[per_index][3]]) {
			for (int seq_index = temp_arr[Constants::threeOptPermutation[per_index][2]];seq_index <= temp_arr[Constants::threeOptPermutation[per_index][3]];seq_index++)
				temp.addNode(seq[seq_index], counter++, g);
		}
		else {
			for (int seq_index = temp_arr[Constants::threeOptPermutation[per_index][2]];seq_index >= temp_arr[Constants::threeOptPermutation[per_index][3]];seq_index--)
				temp.addNode(seq[seq_index], counter++, g);
		}

		for (int seq_index = k; seq_index < V;seq_index++) 
			temp.addNode(seq[seq_index], counter++, g);

		if (per_index == 0 || bestTour.gettotalDistance() > temp.gettotalDistance())
			bestTour.Copy(temp);
	}

	if (this->totalDistance > bestTour.gettotalDistance())
		this->Copy(bestTour);
}

string Tour::printTour() {
	if (lastIndex < 0)
		return "Empty Tour\n";

	string s = "";

	//cout << s << endl;
	if (lastIndex < V - 1)
		s.append("Incompleted Tour:");
	else
		s.append("Completed Tour:");

	//cout << s << endl;
	for (int i = 0; i <= lastIndex;i++) {
		s.append("node ");
		s.append(to_string(seq[i]));
		s.append("-");
		s.append(to_string(int(weight[i])));
		s.append("-");
		//cout << s << endl;
	}
	s.append("node");
	s.append(to_string(seq[0]));
	s.append("; Total Dist = ");
	s.append(to_string(int(totalDistance)));
	s.append("\n");

	return s;
}


void Tour::rvereseArray(int arr[], int start, int end)
{
	while (start < end)
	{
		int temp = arr[start];
		arr[start] = arr[end];
		arr[end] = temp;
		start++;
		end--;
	}
}