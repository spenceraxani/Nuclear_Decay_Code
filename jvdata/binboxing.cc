#include "binboxing.hh"

int main()
{
	cout << endl << "Running bin boxing algorithm..." << endl << endl;

	// import cps data
	ifstream infile1("binned_data_detrended_norm_with_errors.txt");
	double t, cps, delta;
	vector<double> cps_vector, cps_t_vector;
	
	while (infile1 >> t >> cps >> delta)
	{
		cps_t_vector.push_back(t);
		cps_vector.push_back(cps);
	}

	// import long x-ray flux data
	ifstream infile2("XRAY_5m_Long.txt");
	double flux;
	vector<double> longxray_vector, longxray_t_vector;

	while (infile2 >> t >> flux)
	{
		longxray_t_vector.push_back(t);
		longxray_vector.push_back(flux);
	}

	// import short x-ray flux data
	ifstream infile3("XRAY_5m_Short.txt");
	vector<double> shortxray_vector, shortxray_t_vector;

	while (infile3 >> t >> flux)
	{
		shortxray_t_vector.push_back(t);
		shortxray_vector.push_back(flux);
	}

	// find bin averages for various bin widths from 2 hours to 2 days
	for (double binwidth = 2.0/24.0; binwidth <= 48.0/24.0; binwidth += 2.0/24.0)
	{
		cout << binwidth << endl;

		// do the calculations; FIXME do something with the results
		vector<double> cps_averaged = bin_averages(cps_t_vector, cps_vector, binwidth);
		vector<double> long_averaged = bin_averages(longxray_t_vector, longxray_vector, binwidth);
		vector<double> short_averaged = bin_averages(shortxray_t_vector, shortxray_vector, binwidth);
		print_vector(cps_averaged);
	}

	return 0;
}


