#ifdef MILOTTI_MTS
#include <ANN/ANN.h>					// ANN declarations
#endif

// compile with 
// g++ ann_sample.cpp -I/opt/local/include -L/opt/local/lib -lANN -o test

#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation

#include <iostream>						// I/O library
#include <fstream>						// file I/O

using namespace std;					// make std:: accessible

// parameters 

int				k				= 3;			// number of nearest neighbors
int				dim				= 2;			// dimension
double			eps				= 0;			// error bound
int				maxPts			= 1000;			// maximum number of data points

istream*		dataIn			= NULL;			// input for data points
istream*		queryIn			= NULL;			// input for query points

// I/O routines
// else if (!strcmp(directive,"read_query_pts")) {
// 			cin >> arg;							// input file name
// 			readPts(
// 				query_pts,						// point array
// 				query_size,						// number of points
// 				arg,							// file name
// 				QUERY);							// query points
// 			valid_dirty = ANNtrue;				// validation must be redone
// 		}
bool readPt(istream &in, ANNpoint p)			// read point (false on EOF)
{
	for (int i = 0; i < dim; i++) {
		if(!(in >> p[i])) return false;
	}
	return true;
}

void printPt(ostream &out, ANNpoint p)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}

// main program

int main(int argc, char **argv)
{
	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	ifstream dataIn;
	dataIn.open("/home/usersHR/thierry/git_codes/ANN/test/test1-data.pts"); // example from ANN zip

	nPts = 0;									// read data points

	cout << "Data Points:\n";
	while (nPts < maxPts && readPt(dataIn, dataPts[nPts])) {
		printPt(cout, dataPts[nPts]);
		nPts++;
	}

	kdTree = new ANNkd_tree(					// build search structure
					dataPts,					// the data points
					nPts,						// number of points
					dim);						// dimension of space

	double x,y;
	
		cout << "Enter query point: " << endl;	// echo query point
		cout << "x = ";
		cin >> x;
		cout << "y = ";
		cin >> y;
		
		queryPt[0] = x;
		queryPt[1] = y;
		


		kdTree->annkSearch(						// search
				queryPt,						// query point
				k,								// number of near neighbors
				nnIdx,							// nearest neighbors (returned)
				dists,							// distance (returned)
				eps);							// error bound

		cout << "\tNN:\tIndex\tDistance\n";
		for (int i = 0; i < k; i++) {			// print summary
			dists[i] = sqrt(dists[i]);			// unsquare distance
			cout << "\t" << i << "\t" << nnIdx[i] << "\t" << dists[i] << "\n";
			cout << "Point: ";
			printPt(cout, dataPts[nnIdx[i]]);
		}
		
	
    delete [] nnIdx;							// clean things up
    delete [] dists;
    delete kdTree;
	annClose();									// done with ANN

	return EXIT_SUCCESS;
}
