#include "binboxing.hh"

int main()
{
	cout << endl << "Running bin boxing algorithm..." << endl << endl;
	cout << fixed << setprecision(8);

	// import xray spike times
	import_xray_spikes();

	// import cps data
	ifstream infile1("binned_data_detrended_norm_with_errors.txt");
	assert(infile1);
	double t, cps, delta;
	vector<double> cps_vector, cps_t_vector;
	
	while (infile1 >> t >> cps >> delta)
	{
		cps_t_vector.push_back(t);
		cps_vector.push_back(cps);
	}

	double most_unlikely_probability = 1.0; // initialize to most likely and whittle it down
	double corresponding_binwidth = 0.0;
	int corresponding_flare = 0;

	vector<int> zero_prob_flares;
	vector<double> zero_prob_binwidths;
	vector<double> zero_prob_probabilities;

	// loop over the 5 solar flares
	for (int flare_index = 0; flare_index < the_big_five_dates.size(); flare_index++)
	{
		// find bin averages for various bin widths from 2 hours to 2 days
		for (double binwidth = 2.0/24.0; binwidth <= 48.0/24.0; binwidth += 2.0/24.0)
		{
			cout << "solar flare " << flare_index+1 << ", " << "binwidth = " << binwidth << " days" << endl;

			// do the calculations; FIXME do something with the results
			vector<double> cps_averaged = bin_averages(cps_t_vector, cps_vector, binwidth);

			const double threshold =
				average_before_time(cps_t_vector, cps_vector, binwidth, the_big_five_dates[flare_index]);

			cout << " threshold = " << threshold << endl;

			// compute probability that the histogram will contain the threshold value
			double probability;
			if (threshold > 1.00)
				probability = fraction_entries_above(threshold, cps_averaged);
			else if (threshold < 1.00 and threshold > 0.00)
				probability = fraction_entries_below(threshold, cps_averaged);
			else { cerr << "Error with probability calculation." << endl; exit(1); }
			cout << " probability = " << probability << endl;

			// check if we get a new most unlikely probability and modify accordingly
			if (probability < most_unlikely_probability /*and probability > 0.00001*/)
			{
				most_unlikely_probability = probability;
				corresponding_binwidth = binwidth;
				corresponding_flare = flare_index;
			}
			
			if (probability < tol)
			{
				zero_prob_binwidths.push_back(binwidth);
				zero_prob_flares.push_back(flare_index);
				zero_prob_probabilities.push_back(probability);
			}

		} // end binwidth loop
	} // end solar flare loop

	cout << endl;
	cout << "Most unlikely probability = " << most_unlikely_probability << ", corresponding to " <<
		"solar flare " << corresponding_flare+1 << " with binwidth " << corresponding_binwidth <<
		" days" << endl;

	assert(zero_prob_flares.size() == zero_prob_binwidths.size());
	if (zero_prob_flares.size() > 1)
	{
		cout << "More than one zero-probability result (with tolerance " 
				 << tol << ") identified:" << endl;

		cout << "  Solar flares "; print_vector(zero_prob_flares);
		cout << "  Binwidths "; print_vector(zero_prob_binwidths);
		cout << "  Probabilities "; print_vector(zero_prob_probabilities);
	}


	return 0;
}


