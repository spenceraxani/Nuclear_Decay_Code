#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

using namespace std;

// floating-point comparison tolerance
const double tol = 1.0e-8;

// print elements of a vector
template <typename T>
void print_vector(vector<T> vec)
{
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << " ";
	cout << endl;
}

// calculate averages for each bin
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

// five x-ray spikes
vector<double> the_big_five_dates;
void import_xray_spikes()
{
	the_big_five_dates.push_back(56590.0480361);
	the_big_five_dates.push_back(56590.345276);
	the_big_five_dates.push_back(56594.6212844999);
	the_big_five_dates.push_back(56601.6467622);
	the_big_five_dates.push_back(56664.4854034);
}

// calculate the average of yvec within binwidth before final_time
double average_before_time(
	vector<double> tvec, vector<double> yvec, double binwidth, double final_time)
{
	double result = 0.0;
	int counter = 0;
	for (int i = 0; i < tvec.size(); i++)
	{
		const double t = tvec[i];
		if (t >= final_time-binwidth and t < final_time)
		{
			result += yvec[i];
			counter++;
		}
	}
	return (double) result/counter;
}

// fraction of entries above threshold, essentially integrating the distribution
double fraction_entries_above(double threshold, vector<double> yvec)
{
	int result = 0;
	for (int i = 0; i < yvec.size(); i++)
		if (yvec[i] > threshold) result++;
	cout << " yvec.size() = " << yvec.size() << endl;
	return (double) result/yvec.size();
}

// fraction of entries below threshold, essentially integrating the distribution
double fraction_entries_below(double threshold, vector<double> yvec)
{
	int result = 0;
	for (int i = 0; i < yvec.size(); i++)
		if (yvec[i] < threshold) result++;
	cout << " result = " << result << ", yvec.size() = " << yvec.size() << endl;
	return (double) result/yvec.size();
}

