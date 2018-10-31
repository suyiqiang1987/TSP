#include "stdafx.h"
#include "ConstructiveHeuristics.h"

//NNGreedyHeuristic is the nearest neighbourhood greedy heuristic, each iteration add the edge from the previous node with the minimum distance
Tour ConstructiveHeuristics::NNGreedyHeuristic(int startingIndex, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start NN Greedy Heuristic..." << endl;
	begin = clock();

	int V = g.numberOfNodes();
	Tour t(V);
	t.addNode(startingIndex, 0, g);
	for (int i = 1; i < V; i++) {
		//cout << t.printTour();
		int prev_node = t.getSeq(i - 1);
		//cout << "Prev Node" << prev_node << endl;
		double leastEdgeVal = INF;
		double leastEdgeIndex = -1;
		for (int j = 0; j < V;j++) {
			if (t.getHasBeenSeen(j))
				continue;
			if (prev_node == j)
				continue;

			double dist = g.getEdge(prev_node, j);
			//cout << prev_node << " - " << j << ":" << dist <<" v.s. " << leastEdgeVal << endl;;
			if (dist < leastEdgeVal) {
				leastEdgeIndex = j;
				leastEdgeVal = dist;
			}
		}
		t.addNode(leastEdgeIndex, i, g);
	}

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish NN Greedy Heuristic in " << timeSec << "s." << endl;
	return t;
}

//Enumerate all possible permutation of the nodes, and return the cheapest sequence as the resulting tour
Tour ConstructiveHeuristics::Enumeration(int startingIndex, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start Enumeration..." << endl;
	begin = clock();

	int V = g.numberOfNodes();
	int *index_array = new int[V - 1];
	int counter = 0;
	for (int i = 0; i < V;i++) {
		if (startingIndex == i)
			continue;
		index_array[counter] = i;
		counter++;
	}

	Tour t(V);
	sort(index_array, index_array + (V - 1));
	do {
		Tour temp_t(V);
		temp_t.addNode(startingIndex, 0, g);
		for (int i = 0; i < V - 1;i++) {
			temp_t.addNode(index_array[i], i + 1, g);
		}
		if (t.gettotalDistance() == 0 || temp_t.gettotalDistance() < t.gettotalDistance())
		{
			t.Copy(temp_t);
			//cout <<t.printTour();
		}

	} while (next_permutation(index_array, index_array + (V - 1)));

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish Enumeration in " << timeSec << "s." << endl;
	delete[] index_array;
	return t;
}

//Cheapest heuristic insert a node to a sequence which has the lowest cost
Tour ConstructiveHeuristics::CheapestInseration(int startingIndex, Graph& g) {
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start Cheapest Insertion..." << endl;
	begin = clock();

	int V = g.numberOfNodes();
	Tour t(V);
	t.addNode(startingIndex, 0, g);
	for (int iter = 0; iter < V - 1;iter++) {
		double cheapestInsertionVal = INF;
		double cheapestInsertionK = -1;
		double cheapestInserationLocation = -1;
		for (int k = 0; k < V;k++) {
			if (t.getHasBeenSeen(k))
				continue;
			for (int location = 1; location <= t.getLastIndex() + 1;location++) {
				//cout << k << " " << location << endl;
				double temp = t.calculateAddNodeAt(k, location, g);

				if (temp < cheapestInsertionVal) {
					cheapestInsertionVal = temp;
					cheapestInsertionK = k;
					cheapestInserationLocation = location;
				}
			}
		}
		t.addNode(cheapestInsertionK, cheapestInserationLocation, g);
		cout << "..." << t.printTour();
	}

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish Cheapest Inseration in " << timeSec << "s." << endl;
	return t;
}