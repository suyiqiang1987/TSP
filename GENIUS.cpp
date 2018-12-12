#include "stdafx.h"
#include "GENIUS.h"
#include "Constants.h"
#include <algorithm>
#include <iostream>
#include <ctime>

GENIUS::GENIUS(Graph& g, int roriginIndex, int rp, int rtotalNumberOfiterations):bestTour(g.numberOfNodes()){
	clock_t begin;
	clock_t end;
	double timeSec;
	cout << "Start GENIUS..." << endl;
	begin = clock();

	this->V = g.numberOfNodes();
	this->originIndex = roriginIndex;
	this->p = rp;
	this->totalNumberOfiterations = min(rtotalNumberOfiterations, this->V);
	Initialization(g);
	cout << "After initialization:" << bestTour.printTour();

	
		GENI(g);
		//cout << "After GENI:" << bestTour.printTour();

		US(g);
		cout << "Final:" << bestTour.printTour();
	

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish GENIUS in " << timeSec << "s." << endl;

}
void GENIUS::Initialization(Graph& g) {
	auto rng = default_random_engine{};

	this->pNeighbourhood = new int*[this->V];
	this->edges = new double*[this->V];
	for (int i = 0; i < V;i++) {
		this->pNeighbourhood[i] = new int[p];
		this->edges[i] = new double[this->V];
		for (int j = 0; j < p;j++) {
			pNeighbourhood[i][j] = -1;
		}
		for (int j = 0; j < V;j++)
			this->edges[i][j] = g.getEdge(i, j);

		if (i != this->originIndex)
			this->unsolvedNode.push_back(i);
	}

	this->bestTour.addNode(this->originIndex, 0, g);
	updatePNeighbourhood(this->originIndex, g);


	for (int i = 1; i <= min(3,this->V-1); i++) {
		shuffle(begin(unsolvedNode), end(unsolvedNode), rng);
		int temp_index = unsolvedNode.back();
		unsolvedNode.pop_back();
		bestTour.addNode(temp_index, i, g);
		updatePNeighbourhood(temp_index, g);
	}

}

void GENIUS::updatePNeighbourhood(int newIndex, Graph& g) {
	for (int i = 0; i < V;i++) {
		if (i == newIndex)
			continue;

		int largetIndex = 0;
		double largestValue = -1;
		for (int j = 0; j < p;j++) {
			if (pNeighbourhood[i][j] == -1) {
				largetIndex = -1;
				largestValue = -1;
				pNeighbourhood[i][j] = newIndex;
				break;
			}
			if (g.getEdge(i, pNeighbourhood[i][j])>largestValue) {
				largestValue = g.getEdge(i, pNeighbourhood[i][j]);
				largetIndex = j;
			}
		}

		if (largetIndex >= 0) {
			if (g.getEdge(i, newIndex) < largestValue) {
				pNeighbourhood[i][largetIndex] = newIndex;
			}
		}
	}
}



void GENIUS::GENI(Graph &g) {
	auto rng = default_random_engine{};
	for (int iter = 4;iter < V;iter++) {
		shuffle(begin(unsolvedNode), end(unsolvedNode), rng);
		int temp_index = unsolvedNode.back();
		unsolvedNode.pop_back();

		Tour typeIBestTour(this->V);
		typeIInsertion(g, iter, temp_index, typeIBestTour,bestTour);

		Tour typeIIBestTour(this->V);
		typeIIInsertion(g, iter, temp_index, typeIIBestTour, bestTour);
		
		if(typeIBestTour.gettotalDistance() < typeIIBestTour.gettotalDistance())
			bestTour.Copy(typeIBestTour);
		else
			bestTour.Copy(typeIIBestTour);
		//cout << iter << ":" << bestTour.printTour();

		updatePNeighbourhood(temp_index, g);
		//break;
	}
}
//originall array is in the sequence i,i+1,..., j-1, j
//inverse the sequence to j,j-1,...,i+1,i
void GENIUS::typeIInsertion(Graph& g, int iter, int nodeV, Tour& typeIBestTour, Tour& inputTour) {
	
	double bestTourLength = Constants::INF;

	for (int i = 0; i < iter;i++) {
		int indexI = inputTour.getSeq(i);
		int indexIPlusOne = inputTour.getSeq((i+1)%iter);
		for (int j = i + 1;j <= i + iter - 2;j++) {
			int indexJ = inputTour.getSeq(j%iter);
			if (isInPNeighbourHood(nodeV, indexJ) == false) // vj must be pNeighbour of node V
				continue;
			for (int k = j + 1; k <=i + iter - 1;k++) {
				int indexK = inputTour.getSeq(k%iter);
				if (isInPNeighbourHood(indexIPlusOne, indexK) == false)//vk must be pNeighbour of node vi+1
					continue;
				Tour temp(this->V);
				int* temp_array = new int[iter+1];
				int conter = 0;
				temp_array[conter++] = nodeV; //add node (vi,v);
				reverse((i + 1)%iter, j%iter, temp_array, conter, inputTour);//reverse the seqeuqnce (vi+1,...,vj)
				reverse((j + 1)%iter, k%iter, temp_array, conter, inputTour);//reverse the seqeuqnce (vi+1,...,vj)
				addBetween((k + 1)%iter, i, temp_array, conter, inputTour);//add sequence(vk+1,...,vi);
				int temp_array_origin_index = -1;
				for (int l = 0; l < iter + 1;l++) {
					if (temp_array[l] == this->originIndex)
					{
						temp_array_origin_index = l;
					}
				}
				for (int l = 0;l < iter + 1;l++) {
					int indexL = (l + temp_array_origin_index) % (iter + 1);
					temp.addNode(temp_array[indexL], l, g);
				}
				//cout << temp.printTour();

				if (temp.gettotalDistance() < bestTourLength) {
					typeIBestTour.Copy(temp);
					bestTourLength = temp.gettotalDistance();
				}
				
				delete [] temp_array;
			}
		}
	}
}
void GENIUS::typeIIInsertion(Graph& g, int iter, int nodeV, Tour& typeIIBestTour, Tour& inputTour) {

	int* seq = new int[iter + 1];
	int* __restrict cur_seq = new int[iter];
	for (int index = 0; index < iter;index++)
		cur_seq[index] = inputTour.getSeq(index);

	double bestTourLength = Constants::INF;
	//,copyin(cur_seq),copy(seq),copy(temp_array)
   // #pragma acc data copyin(cur_seq),copy(temp_array),copy(seq)
	{
		for (int i = 0; i < iter;i++) {
			//int indexI = inputTour.getSeq(i);
			//int indexI = cur_seq[i];
			//int indexIPlusOne = inputTour.getSeq((i + 1) % iter);
			int indexIPlusOne =cur_seq[(i + 1) % iter];
			for (int l = i + 2; l <= i + iter - 2;l++) {
				//int indexL = inputTour.getSeq(l%iter);
				int indexL = cur_seq[l%iter];
				for (int j = l; j <= i + iter - 2;j++) {
					//int indexJ = inputTour.getSeq(j%iter);
					int indexJ = cur_seq[j%iter];
					//int indexJPlusOne = inputTour.getSeq((j + 1) % iter);
					int indexJPlusOne = cur_seq[(j + 1) % iter];
					if (isInPNeighbourHood(indexJPlusOne, indexL) == false)//vl must be pNeighbour of node vj+1
						continue;
					if (isInPNeighbourHood(nodeV, indexJ) == false) // vj must be pNeighbour of node V
						continue;
					for (int k = j + 2; k <= i + iter;k++) {
						//int indexK = inputTour.getSeq(k%iter);
						int indexK =cur_seq[k%iter];
						if (isInPNeighbourHood(indexIPlusOne, indexK) == false)//vk must be pNeighbour of node vi+1
							continue;

						int conter = 0;
						int* temp_array = new int[iter + 1];
						temp_array[conter++] = nodeV; //add node (vi,v);
						//reverse(l % iter, j%iter, temp_array, conter, inputTour);//reverse the seqeuqnce (vi+1,...,vj)
						reverse(l % iter, j%iter, temp_array, conter, cur_seq,iter-1);//reverse the seqeuqnce (vi+1,...,vj)
						//addBetween((j + 1) % iter, (k - 1) % iter, temp_array, conter, inputTour);//add sequence (vj+1,..., vk-1)
						addBetween((j + 1) % iter, (k - 1) % iter, temp_array, conter, cur_seq, iter - 1);//add sequence (vj+1,..., vk-1)
						//reverse((i + 1) % iter, (l - 1) % iter, temp_array, conter, inputTour);//reverse the seqeuqnce (vi+1,...,vl-1)
						reverse((i + 1) % iter, (l - 1) % iter, temp_array, conter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vl-1)
						//addBetween(k % iter, i, temp_array, conter, inputTour);//add sequence(vk+1,...,vi);
						addBetween(k % iter, i, temp_array, conter, cur_seq, iter - 1);//add sequence(vk+1,...,vi);

						int temp_array_origin_index = -1;
						for (int m = 0; m < iter + 1;m++) {
							if (temp_array[m] == this->originIndex)
							{
								temp_array_origin_index = m;
							}
						}

						double length = 0;
						for (int m = 0; m < iter + 1;m++) {
							int indexM = (m + temp_array_origin_index) % (iter + 1);
							int indexMPlusOne = (m + temp_array_origin_index + 1) % (iter + 1);
							length += this->edges[temp_array[indexM]][temp_array[indexMPlusOne]];
							//temp.addNode(temp_array[indexM], m, g);
						}

						//tours.push_back(temp);
						if (length < bestTourLength) {
							//typeIIBestTour.Copy(temp);
							for (int m = 0; m < iter + 1;m++) {
								int indexM = (m + temp_array_origin_index) % (iter + 1);
								seq[m] = temp_array[indexM];
								//cout << seq[i] << endl;
							}
							bestTourLength = length;
						}

						delete[] temp_array;

					}
				}
			}
		}
	}
	for (int i = 0; i < iter + 1;i++)
		typeIIBestTour.addNode(seq[i], i, g);


	delete[] cur_seq;
	delete[] seq;
	

}


void GENIUS::US(Graph& g) {//Improve Phase of GENIUS
	int t = 1;
	int i = 1;
	double tau_start = bestTour.gettotalDistance();
	while (t <= this->totalNumberOfiterations) {
		Tour iterTour(this->V);
		iterTour.Copy(bestTour);

		Tour typeIUNstringBestTour(this->V);
		typeIUNstringBestTour.Copy(bestTour);
		typeIUnstring(g, i, typeIUNstringBestTour, iterTour);

		Tour typeIIUNstringBestTour(this->V);
		typeIIUNstringBestTour.Copy(bestTour);
		typeIIUnstring(g, i, typeIIUNstringBestTour, iterTour);

		if (typeIUNstringBestTour.gettotalDistance() < typeIIUNstringBestTour.gettotalDistance())
			iterTour.Copy(typeIUNstringBestTour);
		else
			iterTour.Copy(typeIIUNstringBestTour);

		if (iterTour.gettotalDistance() < bestTour.gettotalDistance()) {
			bestTour.Copy(iterTour);
			//cout << t << ":" << bestTour.printTour();
			i = 1;
			t = 1;
		}
		else {
			i++;
			t++;
		}
	}
}


void GENIUS::typeIUnstring(Graph& g, int i, Tour& typeIUNstringBestTour, Tour& t) {
	if (i <=0 || i>=V)//The origin cannt be dropped
	{
		return;
	}


	int IndexI = t.getSeq(i);
	int IndexIPlusOne = t.getSeq((i + 1) % V);
	int indexMinusOne = t.getSeq((i - 1) % V);
	//cout << IndexI << " - " << IndexIPlusOne << " - " << indexMinusOne << endl;
	{
		for (int j = i + 2; j <= i + V - 2;j++) {
			int indexJ = t.getSeq(j%this->V);
			//cout << indexJ << endl;
			if (isInPNeighbourHood(IndexIPlusOne, indexJ) == false) // vj must be pNeighbour of node vi+1
				continue;

			for (int k = i + 1;k <= j - 1;k++) {
				int indexK = t.getSeq(k%this->V);
				//cout << indexK << endl;
				if (isInPNeighbourHood(indexMinusOne, indexK) == false) // vk must be pNeighbour of node vk-1
					continue;

				Tour temp(this->V);
				int* temp_array = new int[this->V - 1];
				int conter = 0;

				//cout << i << " " << j << " " << k << endl;;
				reverse((i + 1) % this->V, k%this->V, temp_array, conter, t);//reverse the seqeuqnce (vi+1,...,vj)
				reverse((k + 1) % this->V, j%this->V, temp_array, conter, t);//reverse the seqeuqnce (vi+1,...,vj)
				addBetween((j + 1) % this->V, (i - 1) % this->V, temp_array, conter, t);//add sequence(vk+1,...,vi);
				int temp_array_origin_index = -1;
				for (int l = 0; l < this->V - 1;l++) {
					if (temp_array[l] == this->originIndex)
					{
						temp_array_origin_index = l;
					}
				}
				for (int l = 0;l < this->V - 1;l++) {
					int indexL = (l + temp_array_origin_index) % (this->V - 1);
					temp.addNode(temp_array[indexL], l, g);
				}

				//cout << "Temp Tour:" << temp.printTour() << endl;

				Tour typeIBestTour(this->V);
				typeIInsertion(g, this->V - 1, IndexI, typeIBestTour, temp);
				//cout << "Type I Insertion Tour:" << typeIBestTour.printTour() << endl;

				Tour typeIIBestTour(this->V);
				typeIIInsertion(g, this->V - 1, IndexI, typeIIBestTour, temp);
				//cout << "Type II Insertion Tour:" << typeIIBestTour.printTour() << endl;

				if (typeIBestTour.gettotalDistance() < typeIIBestTour.gettotalDistance())
					temp.Copy(typeIBestTour);
				else
					temp.Copy(typeIIBestTour);

				if (temp.gettotalDistance() < typeIUNstringBestTour.gettotalDistance())
					typeIUNstringBestTour.Copy(temp);

				//cout << "here" << endl;
				delete[] temp_array;
				//cout << "there" << endl;


			}
		}
	}
}

void GENIUS::typeIIUnstring(Graph& g, int i, Tour& typeIIUNstringBestTour, Tour& t) {
	if (i <= 0 || i >= V)//The origin cannt be dropped
	{
		return;
	}

	int IndexI = t.getSeq(i);
	int IndexIPlusOne = t.getSeq((i + 1) % V);
	int indexMinusOne = t.getSeq((i - 1) % V);

	for (int j = i + 2; j <= i + V - 3;j++) {
		int indexJ = t.getSeq(j%this->V);
		//cout << indexJ << endl;
		if (isInPNeighbourHood(IndexIPlusOne, indexJ) == false) // vj must be pNeighbour of node vi+1
			continue;

		for (int l = j; l <= i+V - 3;l++) {
			int indexL = t.getSeq(l%this->V);

			for (int k = l +1;k <= i+V-2;k++) {
				int indexK = t.getSeq(k%this->V);
				int indexKPlusOne = t.getSeq((k+1)%this->V);
				//cout << indexK << endl;
				if (isInPNeighbourHood(indexMinusOne, indexK) == false) // vk must be pNeighbour of node vk-1
					continue;

				if (isInPNeighbourHood(indexKPlusOne, indexL) == false) // vk must be pNeighbour of node vk-1
					continue;

				Tour temp(this->V);
				int* temp_array = new int[this->V - 1];
				int conter = 0;

				
				reverse((i + 1) % this->V, (j-1)%this->V, temp_array, conter, t);//reverse the seqeuqnce (vi+1,...,vj-1)
				addBetween(j % this->V, l% this->V, temp_array, conter, t);//add sequence(vj,...,vl);
				addBetween((k + 1) % this->V, (i - 1) % this->V, temp_array, conter, t);//add sequence(vk+1,...,vi-1);
				reverse((l + 1) % this->V, k%this->V, temp_array, conter, t);//reverse the seqeuqnce (vl+1,...,vk)
				
				int temp_array_origin_index = -1;
				for (int m = 0; m < this->V - 1;m++) {
					if (temp_array[m] == this->originIndex)
					{
						temp_array_origin_index = m;
					}
				}
				for (int m = 0;m < this->V - 1;m++) {
					int indexM = (m + temp_array_origin_index) % (this->V - 1);
					temp.addNode(temp_array[indexM], m, g);
				}

				Tour typeIBestTour(this->V);
				typeIInsertion(g, this->V - 1, IndexI, typeIBestTour, temp);
				//cout << "Type I Insertion Tour:" << typeIBestTour.printTour() << endl;

				Tour typeIIBestTour(this->V);
				typeIIInsertion(g, this->V - 1, IndexI, typeIIBestTour, temp);
				//cout << "Type II Insertion Tour:" << typeIIBestTour.printTour() << endl;

				if (typeIBestTour.gettotalDistance() < typeIIBestTour.gettotalDistance())
					temp.Copy(typeIBestTour);
				else
					temp.Copy(typeIIBestTour);

				if (temp.gettotalDistance() < typeIIUNstringBestTour.gettotalDistance())
					typeIIUNstringBestTour.Copy(temp);

				//cout << "here" << endl;
				delete[] temp_array;
				//cout << "there" << endl;
			}

		}
	}
}

bool GENIUS::isInPNeighbourHood(int i, int j){
	for (int k = 0; k < p;k++) {
		if (pNeighbourhood[i][k] == j)
			return true;
	}
	return false;
}

void GENIUS::reverse(int i, int j, int* temp_array, int& counter, Tour t) {
	if (i <= j) {
		for (int k = j; k >= i;k--) {
			temp_array[counter++] = t.getSeq(k);
		}
	}
	else {
		for (int k = j; k >= 0;k--) {
			temp_array[counter++] = t.getSeq(k);
		}
		for (int k = t.getLastIndex(); k >= i;k--) {
			temp_array[counter++] = t.getSeq(k);
		}
	}
}


void GENIUS::reverse(int i, int j, int* temp_array, int& counter, int* t, int lastIndex) {
	if (i <= j) {
		for (int k = j; k >= i;k--) {
			temp_array[counter++] = t[k];
		}
	}
	else {
		for (int k = j; k >= 0;k--) {
			temp_array[counter++] = t[k];
		}
		for (int k = lastIndex; k >= i;k--) {
			temp_array[counter++] = t[k];
		}
	}
}

void GENIUS::addBetween(int i, int j, int* temp_array, int& counter, Tour t) {
	if (i <= j) {
		for (int k = i; k <= j;k++) {
			temp_array[counter++] = t.getSeq(k);
		}
	}
	else {
		for (int k = i; k <= t.getLastIndex();k++) {
			temp_array[counter++] = t.getSeq(k);
		}
		for (int k = 0; k <= j;k++) {
			temp_array[counter++] = t.getSeq(k);
		}
	}
}

void GENIUS::addBetween(int i, int j, int* temp_array, int& counter, int *t, int lastIndex) {
	if (i <= j) {
		for (int k = i; k <= j;k++) {
			temp_array[counter++] = t[k];
		}
	}
	else {
		for (int k = i; k <= lastIndex;k++) {
			temp_array[counter++] = t[k];
		}
		for (int k = 0; k <= j;k++) {
			temp_array[counter++] = t[k];
		}
	}
}

