#include "stdafx.h"
#include "GENIUS_GPU.h"
#include "Constants.h"
#include <algorithm>
#include <iostream>
#include <ctime>
const int GPUThreshold = 32;
GENIUS_GPU::GENIUS_GPU(Graph& g, int roriginIndex, int rp, int rtotalNumberOfiterations):bestTour(g.numberOfNodes()){
	clock_t begin;
	clock_t end;
	double timeSec;
	this->V = g.numberOfNodes();
	double *__restrict*__restrict edges = new double*__restrict[this->V];
	for (int i = 0; i < V;i++) {
		edges[i] = new double[this->V];
		for (int j = 0; j < V;j++)
			edges[i][j] = g.getEdge(i, j);
	}


	cout << "Start GENIUS..." << endl;
	begin = clock();


	this->originIndex = roriginIndex;
	this->p = rp;
	this->totalNumberOfiterations = min(rtotalNumberOfiterations, this->V);
	Initialization(g);
	cout << "After initialization:" << bestTour.printTour();

	{
		GENI(g,edges, pNeighbourhood);
		//cout << "After GENI:" << bestTour.printTour();

		US(g, edges, pNeighbourhood);
		cout << "Final:" << bestTour.printTour();
	}

	end = clock();
	timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "...Finish GENIUS in " << timeSec << "s." << endl;

	delete[] edges;

}
void GENIUS_GPU::Initialization(Graph& g) {
	auto rng = default_random_engine{};

	this->pNeighbourhood = new int*[this->V];
	for (int i = 0; i < V;i++) {
		this->pNeighbourhood[i] = new int[p];
		for (int j = 0; j < p;j++) {
			pNeighbourhood[i][j] = -1;
		}
		
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

void GENIUS_GPU::updatePNeighbourhood(int newIndex, Graph& g) {
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



void GENIUS_GPU::GENI(Graph &g, double *__restrict*__restrict e, int**  pNeighbourhood) {
	auto rng = default_random_engine{};
	for (int iter = 4;iter < V;iter++) {
		shuffle(begin(unsolvedNode), end(unsolvedNode), rng);
		int temp_index = unsolvedNode.back();
		unsolvedNode.pop_back();

		Tour typeIBestTour(this->V);
		if(iter<= GPUThreshold)
			typeIInsertion(g, iter, temp_index, typeIBestTour,bestTour,e, pNeighbourhood);
		else
			typeIInsertionGPU(g, iter, temp_index, typeIBestTour, bestTour, e, pNeighbourhood);

		Tour typeIIBestTour(this->V);
		typeIIInsertion(g, iter, temp_index, typeIIBestTour, bestTour,e, pNeighbourhood);
		
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
void GENIUS_GPU::typeIInsertion(Graph& g, const int iter, int nodeV, Tour& typeIBestTour, Tour& inputTour, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {
	int* __restrict cur_seq = new int[iter];
	for (int index = 0; index < iter;index++)
		cur_seq[index] = inputTour.getSeq(index);
	double bestTourLength = Constants::INF;
	int bestI = -1;
	int bestJ = -1;
	int bestK = -1;

	for (int i = 0; i < iter;i++) {
		int indexIPlusOne = cur_seq[(i + 1) % iter];
   
		for (int j = i+1; j <= i+iter - 2;j++){
				
			int indexJ = cur_seq[j%iter];
			if (isInPNeighbourHood(nodeV, indexJ,p, pNeighbourhood) == false) // vj must be pNeighbour of node V
				continue;
               
			for (int k = j+1; k <= i+iter - 1;k++){

				int indexK = cur_seq[k%iter];
				if (isInPNeighbourHood(indexIPlusOne, indexK,p, pNeighbourhood) == false)//vk must be pNeighbour of node vi+1
					continue;

				double length = 0;
				length += e[nodeV][(cur_seq[j%iter])];
				length += reverseDist((i + 1) % iter, j%iter,  cur_seq, iter - 1,e);
				length += e[(cur_seq[(i+1)%iter])][(cur_seq[(k) % iter])];
				length += reverseDist((j + 1) % iter, k%iter, cur_seq, iter - 1,e);
				length += e[(cur_seq[(j + 1) % iter])][(cur_seq[(k+1) % iter])];
				length += addBetweenDist((k + 1) % iter, i, cur_seq, iter - 1,e);
				length += e[(cur_seq[(i) % iter])][nodeV];
				if (length < bestTourLength) {
					bestI = i;
					bestJ = j;
					bestK = k;
					bestTourLength = length;
				}
			}
		}
	}

	//cout << bestI << bestJ << bestK << endl;
	int counter = 0;
	int* seq = new int[iter + 1];
	seq[counter++] = nodeV;
	reverse((bestI + 1) % iter, bestJ%iter, seq, counter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vj)
	reverse((bestJ + 1) % iter, bestK%iter, seq, counter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vj)
	addBetween((bestK + 1) % iter, bestI, seq, counter, cur_seq, iter - 1);//add sequence(vk+1,...,vi);
	int temp_array_origin_index = -1;
	for (int m = 0; m < iter + 1;m++) {
		if (seq[m] == this->originIndex)
		{
			temp_array_origin_index = m;
		}
		//cout << m << ":" << seq[m] << endl;
	}
	Tour temp(this->V);
	for (int m = 0; m < iter + 1;m++) {
		int indexM = (m + temp_array_origin_index) % (iter + 1);
		temp.addNode(seq[indexM], m, g);
	}
	typeIBestTour.Copy(temp);
	delete[]seq;
	delete[] cur_seq;
	
}


void GENIUS_GPU::typeIInsertionGPU(Graph& g, const int iter, int nodeV, Tour& typeIBestTour, Tour& inputTour, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {
	int* __restrict cur_seq = new int[iter];
	for (int index = 0; index < iter;index++)
		cur_seq[index] = inputTour.getSeq(index);
	double bestTourLength = Constants::INF;

	double *__restrict dist_matrix = new double[V*V*V];
	//#pragma acc region
    #pragma acc parallel loop device_type(nvidia) gang worker num_workers(V)
	{
		#pragma acc loop independent 
		for (int i = 0; i < V;i++)
		{
			#pragma acc loop independent
			for (int j = 0; j < V;j++) {
                //#pragma acc loop independent
				for (int k = 0; k < V;k++) {
					dist_matrix[i*V*V + j*V + k] = Constants::INF;
				}
			}
		}
	}
    #pragma acc data copy(dist_matrix[0:V][0:V][0:V]), copyin(e[0:V][0:V],cur_seq[0:iter],this,iter,nodeV,p)
	{
		//#pragma acc for private(dist_matrix[0:50][0:50][0:50])
        #pragma acc parallel loop device_type(nvidia) gang worker num_workers(iter)
		for (int i = 0; i < iter;i++) {
			int indexIPlusOne = cur_seq[(i + 1) % iter];

            #pragma acc loop independent
			//for (int j = i + 1;j <= i + iter - 2;j++) {
			for (int j = i + 1; j <= i + iter - 2;j++) {
				//cout << j << endl;
				int indexJ = cur_seq[j%iter];
				if (isInPNeighbourHood(nodeV, indexJ, p, pNeighbourhood) == false) // vj must be pNeighbour of node V
					continue;

				//for (int k = j + 1; k <= i + iter - 1;k++) {
                #pragma acc loop independent
				for (int k = j + 1; k <= i + iter - 1;k++) {

					int indexK = cur_seq[k%iter];
					if (isInPNeighbourHood(indexIPlusOne, indexK, p, pNeighbourhood) == false)//vk must be pNeighbour of node vi+1
						continue;


					double length = 0;
					length += e[nodeV][(cur_seq[j%iter])];
					length += reverseDist((i + 1) % iter, j%iter, cur_seq, iter - 1, e);
					length += e[(cur_seq[(i + 1) % iter])][(cur_seq[(k) % iter])];
					length += reverseDist((j + 1) % iter, k%iter, cur_seq, iter - 1, e);
					length += e[(cur_seq[(j + 1) % iter])][(cur_seq[(k + 1) % iter])];
					length += addBetweenDist((k + 1) % iter, i, cur_seq, iter - 1, e);
					length += e[(cur_seq[(i) % iter])][nodeV];

					dist_matrix[(i%iter)*V*V + (j%iter)*V + k%iter] = length;

				}
			}
		}
	}


	int bestI = -1;
	int bestJ = -1;
	int bestK = -1;
	//#pragma acc parallel loop
	{
		for (int i = 0; i < V;i++) {
			for (int j = 0; j < V;j++) {
				for (int k = 0; k < V;k++) {
					if (bestTourLength > dist_matrix[i*V*V + j*V + k]) {
						bestTourLength = dist_matrix[i*V*V + j*V + k];
						bestI = i;
						bestJ = j;
						bestK = k;
					}
				}
			}
		}
	}
	//cout << bestI << bestJ << bestK << endl;
	int counter = 0;
	int* seq = new int[iter + 1];
	seq[counter++] = nodeV;
	reverse((bestI + 1) % iter, bestJ%iter, seq, counter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vj)
	reverse((bestJ + 1) % iter, bestK%iter, seq, counter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vj)
	addBetween((bestK + 1) % iter, bestI, seq, counter, cur_seq, iter - 1);//add sequence(vk+1,...,vi);
	int temp_array_origin_index = -1;
	for (int m = 0; m < iter + 1;m++) {
		if (seq[m] == this->originIndex)
		{
			temp_array_origin_index = m;
		}
		//cout << m << ":" << seq[m] << endl;
	}
	Tour temp(this->V);
	for (int m = 0; m < iter + 1;m++) {
		int indexM = (m + temp_array_origin_index) % (iter + 1);
		temp.addNode(seq[indexM], m, g);
	}
	typeIBestTour.Copy(temp);
	delete[]seq;
	delete[] cur_seq;

}
void GENIUS_GPU::typeIIInsertion(Graph& g, int iter, int nodeV, Tour& typeIIBestTour, Tour& inputTour, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {

	int* seq = new int[iter + 1];
	int* cur_seq = new int[iter];
	for (int index = 0; index < iter;index++)
		cur_seq[index] = inputTour.getSeq(index);

	double bestTourLength = Constants::INF;
	//,copyin(cur_seq),copy(seq),copy(temp_array)
   // #pragma acc data copyin(cur_seq),copy(temp_array),copy(seq)
   // #pragma acc data present(edges)
    //#pragma acc kernels
	{
		for (int i = 0; i < iter;i++) {
			int indexIPlusOne =cur_seq[(i + 1) % iter];
			for (int l = i + 2; l <= i + iter - 2;l++) {
				int indexL = cur_seq[l%iter];
				for (int j = l; j <= i + iter - 2;j++) {
					int indexJ = cur_seq[j%iter];
					int indexJPlusOne = cur_seq[(j + 1) % iter];
					if (isInPNeighbourHood(indexJPlusOne, indexL,p, pNeighbourhood) == false)//vl must be pNeighbour of node vj+1
						continue;
					if (isInPNeighbourHood(nodeV, indexJ, p,pNeighbourhood) == false) // vj must be pNeighbour of node V
						continue;
					for (int k = j + 2; k <= i + iter;k++) {
						int indexK =cur_seq[k%iter];
						if (isInPNeighbourHood(indexIPlusOne, indexK,p, pNeighbourhood) == false)//vk must be pNeighbour of node vi+1
							continue;

						int conter = 0;
						int* temp_array = new int[iter + 1];
						temp_array[conter++] = nodeV; //add node (vi,v);
						reverse(l % iter, j%iter, temp_array, conter, cur_seq,iter-1);//reverse the seqeuqnce (vi+1,...,vj)
						addBetween((j + 1) % iter, (k - 1) % iter, temp_array, conter, cur_seq, iter - 1);//add sequence (vj+1,..., vk-1)
						reverse((i + 1) % iter, (l - 1) % iter, temp_array, conter, cur_seq, iter - 1);//reverse the seqeuqnce (vi+1,...,vl-1)
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
							length += e[temp_array[indexM]][temp_array[indexMPlusOne]];
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


void GENIUS_GPU::US(Graph& g, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {//Improve Phase of GENIUS
	int t = 1;
	int i = 1;
	double tau_start = bestTour.gettotalDistance();
	while (t <= this->totalNumberOfiterations) {
		Tour iterTour(this->V);
		iterTour.Copy(bestTour);

		Tour typeIUNstringBestTour(this->V);
		typeIUNstringBestTour.Copy(bestTour);
		typeIUnstring(g, i, typeIUNstringBestTour, iterTour,e, pNeighbourhood);

		Tour typeIIUNstringBestTour(this->V);
		typeIIUNstringBestTour.Copy(bestTour);
		typeIIUnstring(g, i, typeIIUNstringBestTour, iterTour,e, pNeighbourhood);

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


void GENIUS_GPU::typeIUnstring(Graph& g, int i, Tour& typeIUNstringBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {
	if (i <=0 || i>=V)//The origin cannt be dropped
	{
		return;
	}

	int *cur_seq = new int[V];
	for (int index = 0; index < V;index++)
		cur_seq[index] = t.getSeq(index);

	int IndexI = t.getSeq(i);
	int IndexIPlusOne = t.getSeq((i + 1) % V);
	int indexMinusOne = t.getSeq((i - 1) % V);
	//cout << IndexI << " - " << IndexIPlusOne << " - " << indexMinusOne << endl;
	{
		for (int j = i + 2; j <= i + V - 2;j++) {
			int indexJ = t.getSeq(j%this->V);
			//cout << indexJ << endl;
			if (isInPNeighbourHood(IndexIPlusOne, indexJ, p,pNeighbourhood) == false) // vj must be pNeighbour of node vi+1
				continue;

			for (int k = i + 1;k <= j - 1;k++) {
				int indexK = t.getSeq(k%this->V);
				//cout << indexK << endl;
				if (isInPNeighbourHood(indexMinusOne, indexK,p, pNeighbourhood) == false) // vk must be pNeighbour of node vk-1
					continue;

				Tour temp(this->V);
				int* temp_array = new int[this->V - 1];
				int conter = 0;

				//cout << i << " " << j << " " << k << endl;;
				reverse((i + 1) % this->V, k%this->V, temp_array, conter, cur_seq,V-1);//reverse the seqeuqnce (vi+1,...,vj)
				reverse((k + 1) % this->V, j%this->V, temp_array, conter, cur_seq,V-1);//reverse the seqeuqnce (vi+1,...,vj)
				addBetween((j + 1) % this->V, (i - 1) % this->V, temp_array, conter, cur_seq,V-1);//add sequence(vk+1,...,vi);
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
				if(this->V<GPUThreshold)
					typeIInsertion(g, this->V - 1, IndexI, typeIBestTour, temp,e, pNeighbourhood);
				else
					typeIInsertionGPU(g, this->V - 1, IndexI, typeIBestTour, temp, e, pNeighbourhood);
				//cout << "Type I Insertion Tour:" << typeIBestTour.printTour() << endl;

				Tour typeIIBestTour(this->V);
				typeIIInsertion(g, this->V - 1, IndexI, typeIIBestTour, temp,e, pNeighbourhood);
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

	delete[] cur_seq;
}

void GENIUS_GPU::typeIIUnstring(Graph& g, int i, Tour& typeIIUNstringBestTour, Tour& t, double *__restrict*__restrict e, int** __restrict pNeighbourhood) {
	if (i <= 0 || i >= V)//The origin cannt be dropped
	{
		return;
	}

	int *cur_seq = new int[V];
	for (int index = 0; index < V;index++)
		cur_seq[index] = t.getSeq(index);

	int IndexI = t.getSeq(i);
	int IndexIPlusOne = t.getSeq((i + 1) % V);
	int indexMinusOne = t.getSeq((i - 1) % V);

	for (int j = i + 2; j <= i + V - 3;j++) {
		int indexJ = cur_seq[j%this->V];
		//cout << indexJ << endl;
		if (isInPNeighbourHood(IndexIPlusOne, indexJ,p, pNeighbourhood) == false) // vj must be pNeighbour of node vi+1
			continue;

		for (int l = j; l <= i+V - 3;l++) {
			int indexL = cur_seq[l%this->V];

			for (int k = l +1;k <= i+V-2;k++) {
				int indexK = cur_seq[k%this->V];
				int indexKPlusOne = cur_seq[(k+1)%this->V];
				//cout << indexK << endl;
				if (isInPNeighbourHood(indexMinusOne, indexK,p, pNeighbourhood) == false) // vk must be pNeighbour of node vk-1
					continue;

				if (isInPNeighbourHood(indexKPlusOne, indexL,p, pNeighbourhood) == false) // vk must be pNeighbour of node vk-1
					continue;

				Tour temp(this->V);
				int* temp_array = new int[this->V - 1];
				int conter = 0;

				
				reverse((i + 1) % this->V, (j-1)%this->V, temp_array, conter, cur_seq,V-1);//reverse the seqeuqnce (vi+1,...,vj-1)
				addBetween(j % this->V, l% this->V, temp_array, conter, cur_seq,V-1);//add sequence(vj,...,vl);
				addBetween((k + 1) % this->V, (i - 1) % this->V, temp_array, conter, cur_seq,V-1);//add sequence(vk+1,...,vi-1);
				reverse((l + 1) % this->V, k%this->V, temp_array, conter, cur_seq,V-1);//reverse the seqeuqnce (vl+1,...,vk)
				
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
				if (this->V<GPUThreshold)
				    typeIInsertion(g, this->V - 1, IndexI, typeIBestTour, temp,e, pNeighbourhood);
				else
					typeIInsertionGPU(g, this->V - 1, IndexI, typeIBestTour, temp, e, pNeighbourhood);
				//cout << "Type I Insertion Tour:" << typeIBestTour.printTour() << endl;

				Tour typeIIBestTour(this->V);
				typeIIInsertion(g, this->V - 1, IndexI, typeIIBestTour, temp,e, pNeighbourhood);
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

	delete[] cur_seq;
}

#pragma acc routine seq
bool GENIUS_GPU::isInPNeighbourHood(int i, int j, int neighboursize, int** __restrict pNeighbourhood){
	//int* point;
	//point = find(pNeighbourhood[i], pNeighbourhood[i] + neighboursize, j);
	//if (point != pNeighbourhood[i] + neighboursize)
		//return true;
	//else
	for (int k = 0; k < neighboursize;k++) {
		if (pNeighbourhood[i][k] == j)
			return true;
	}
		
	return false;

}

void GENIUS_GPU::reverse(int i, int j, int* temp_array, int& counter, Tour t) {
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


void GENIUS_GPU::reverse(int i, int j, int* temp_array, int& counter, int* t, int lastIndex) {
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

#pragma acc routine seq
double GENIUS_GPU::reverseDist(int i, int j, int *__restrict t, int lastIndex, double *__restrict*__restrict e) {
	double length = 0;
	if (i < j) {
        #pragma acc loop reduction(+:length)
		for (int k = j; k >= i + 1;k--)
			length += e[t[k]][t[k - 1]];
	}
	else if (i>j) {
        #pragma acc loop reduction(+:length)
		for (int k = j; k > 0;k--) {
			length += e[t[k]][t[k - 1]];
		}
		length += e[t[0]][t[lastIndex]];
        #pragma acc loop reduction(+:length)
		for (int k = lastIndex; k >= i+1;k--) {
			length += e[t[k]][t[k - 1]];
		}
	}
	return length;
}

void GENIUS_GPU::addBetween(int i, int j, int* temp_array, int& counter, Tour t) {
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

void GENIUS_GPU::addBetween(int i, int j, int* temp_array, int& counter, int *t, int lastIndex) {
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

#pragma acc routine seq
double GENIUS_GPU::addBetweenDist(int i, int j,  int *__restrict t, int lastIndex, double *__restrict*__restrict e) {
	double length = 0;
	if (i < j) {
        #pragma acc loop reduction(+:length)
		for (int k = i; k < j;k++) {
			length += e[t[k]][t[k + 1]];
		}
	}
	else {
       #pragma acc loop reduction(+:length)
		for (int k = i; k < lastIndex;k++) {
			length += e[t[k]][t[k + 1]];
		}
		length += e[t[lastIndex]][t[0]];
        #pragma acc loop reduction(+:length)
		for (int k = 0; k < j;k++) {
			length += e[t[k]][t[k + 1]];
		}
	}

	return length;
}

