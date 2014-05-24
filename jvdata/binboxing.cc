#include "binboxing.hh"

int main()
{
	cout << "Running bin boxing algorithm..." << endl;

	/*ifstream infile("binned_data_detrended_norm_with_errors.txt");
	double t, cps, delta;

	while (infile >> t >> cps >> delta)
	{
		cps_t_vector.push_back(t);
		cps_vector.push_back(cps);
		cout << cps << endl;
		cout << cps_vector.size() << endl;
	}*/

	vector<double> myvec;
	vector<double> xvec;
	for (int i = 0; i < 10; i++)
	{
		xvec.push_back(i);
		myvec.push_back(i*i);
	}

	vector<double> newvec = bin_averages(xvec,myvec,2);

	for (int j = 0; j < newvec.size(); j++)
	{
		cout << newvec[j] << endl;
	}

	return 0;
}


