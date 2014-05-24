#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

// print elements of a vector
void print_vector(vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << " ";
	cout << endl;
}

// calculate averages
vector<double> bin_averages(vector<double> xvec, vector<double> yvec, double binwidth)
{
	vector<double> outvec;
	vector<int> ptsperbin;

	const int length = xvec.size();
	const double min = xvec[0];
	const double max = xvec[length-1];
	const double width = max-min;
	const int nbins = width/binwidth + 1; // automatically gives int floor but cuts off tail end

	// zero everything just in case
	for (int i = 0; i < nbins; i++)
	{
		outvec.push_back(0.0);
		ptsperbin.push_back(0);
	}

	// fill new vectors
	for (int i = 0; i < length; i++)
	{
		outvec[(xvec[i]-min)/binwidth] += yvec[i];
		ptsperbin[(xvec[i]-min)/binwidth]++;
	}

	// average over points per bin
	for (int i = 0; i < outvec.size(); i++)
		outvec[i] /= (1.0*ptsperbin[i]);
	
	return outvec;
}

