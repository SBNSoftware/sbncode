/**
 * \file TruncMean.h
 *
 * \ingroup 3DMichel
 * 
 * \brief Class def header for a class TruncMean
 *
 * @author david caratelli [davidc@fnal.gov]
 * Written 08/02/2018.
 */

#ifndef TRUNCMEAN_H
#define TRUNCMEAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>

/**
   \class TruncMean
   The truncated mean class allows to compute the following quantities
   1) the truncated mean profile of an ordered vector of values, such as 
   the charge profile along a particle's track.
   To create such a profile use the function CalcTruncMeanProfile()
   2) Get the truncated mean value of a distribution. This function
   iteratively hones in on the truncated mean of a distribution by
   updating the mean and cutting the tails after each iteration.
   For this functionality use CalcIterativeTruncMean()
   doxygen documentation!
*/

static const double kINVALID_FLOAT = std::numeric_limits<double>::max();

class TruncMean{

 public:

  /// Default constructor
  TruncMean(){}

  /// Default destructor
  ~TruncMean(){}

  /**
     @brief Given residual range and dq vectors return truncated local dq.
     Input vectors are assumed to be match pair-wise (nth entry in rr_v
     corresponds to nth entry in dq_v vector).
     Input rr_v values are also assumed to be ordered: monotonically increasing
     or decreasing.
     For every dq value a truncated linear dq value is calculated as follows:
     0) all dq values within a rr range set by the class variable _rad are selected.
     1) the median and rms of these values is calculated.
     2) the subset of local dq values within the range [median-rms, median+rms] is selected.
     3) the resulting local truncated dq is the average of this truncated subset.
     @input std::vector<double> rr_v -> vector of x-axis coordinates (i.e. position for track profile)
     @input std::vector<double> dq_v -> vector of measured values for which truncated profile is requested
     (i.e. charge profile of a track)
     @input std::vector<double> dq_trunc_v -> passed by reference -> output stored here
     @input double nsigma -> optional parameter, number of sigma to keep around RMS for TM calculation
  */
  void CalcTruncMeanProfile(const std::vector<double>& rr_v, const std::vector<double>& dq_v,
			    std::vector<double>& dq_trunc_v, const double& nsigma = 1);

  /**
     @brief Iteratively calculate the truncated mean of a distribution
     @input std::vector<double> v -> vector of values for which truncated mean is asked
     @input size_t nmin -> minimum number of iterations to converge on truncated mean
     @input size_t nmax -> maximum number of iterations to converge on truncated mean
     @input size_t lmin -> minimum number of entries in vector before exiting and returning current value
     @input size_t currentiteration -> current iteration
     @input double convergencelimit -> fractional difference between successive iterations
     under which the iteration is completed, provided nmin iterations have occurred.
     @input nsigma -> number of sigma around the median value to keep when the distribution is trimmed.
   */
  double CalcIterativeTruncMean(std::vector<double> v, const size_t& nmin,
			       const size_t& nmax, const size_t& currentiteration,
			       const size_t& lmin,
			       const double& convergencelimit,
			       const double& nsigma, const double& oldmed = kINVALID_FLOAT);

  /**
     @brief Set the smearing radius over which to take hits for truncated mean computaton.
   */
  void setRadius(const double& rad) { _rad = rad; }

 private:

  double Mean  (const std::vector<double>& v);
  double Median(const std::vector<double>& v);
  double RMS   (const std::vector<double>& v);

  /**
     Smearing radius over which charge from neighboring hits is scanned to calculate local
     truncated mean
   */
  double _rad;

};

double TruncMean::CalcIterativeTruncMean(std::vector<double> v, const size_t& nmin,
		const size_t& nmax, const size_t& currentiteration,
		const size_t& lmin,
		const double& convergencelimit,
		const double& nsigma, const double& oldmed)
{

	auto const& mean = Mean(v);
	auto const& med  = Median(v);
	auto const& rms  = RMS(v);

	// if the vector length is below the lower limit -> return
	if (v.size() < lmin)
		return mean;

	// if we have passed the maximum number of iterations -> return
	if (currentiteration >= nmax)
		return mean;

	// if we passed the minimum number of iterations and the mean is close enough to the old value
	double fracdiff = fabs(med-oldmed) / oldmed;
	if ( (currentiteration >= nmin) && (fracdiff < convergencelimit) )
		return mean;

	// if reached here it means we have to go on for another iteration

	// cutoff tails of distribution surrounding the mean
	// use erase-remove : https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
	// https://stackoverflow.com/questions/17270837/stdvector-removing-elements-which-fulfill-some-conditions
	v.erase( std::remove_if( v.begin(), v.end(), 
				[med,nsigma,rms](const double& x) { return ( (x < (med-nsigma*rms)) || (x > (med+nsigma*rms)) ); }), // lamdda condition for events to be removed
			v.end());

	return CalcIterativeTruncMean(v, nmin, nmax, lmin, currentiteration+1, convergencelimit, nsigma, med);
}

void TruncMean::CalcTruncMeanProfile(const std::vector<double>& rr_v, const std::vector<double>& dq_v,
		std::vector<double>& dq_trunc_v, const double& nsigma)
{

	// how many points to sample 
	int Nneighbor = (int)(_rad * 3 * 2);

	dq_trunc_v.clear();
	dq_trunc_v.reserve( rr_v.size() );

	int Nmax = dq_v.size()-1;

	for (size_t n=0; n < dq_v.size(); n++) {

		// current residual range
		double rr = rr_v.at(n);

		int nmin = n - Nneighbor;
		int nmax = n + Nneighbor;

		if (nmin < 0) nmin = 0;
		if (nmax > Nmax) nmax = Nmax;

		// vector for local dq values
		std::vector<double> dq_local_v;

		for (int i=nmin; i < nmax; i++) {

			double dr = rr - rr_v[i];
			if (dr < 0) dr *= -1;

			if (dr > _rad) continue;

			dq_local_v.push_back( dq_v[i] );

		}// for all ticks we want to scan

		if (dq_local_v.size() == 0 ) {
			dq_trunc_v.push_back( dq_v.at(n) );
			continue;
		}

		// calculate median and rms
		double median = Median(dq_local_v);
		double rms    = RMS(dq_local_v);

		double truncated_dq = 0.;
		int npts = 0;
		for (auto const& dq : dq_local_v) {
			if ( ( dq < (median+rms * nsigma) ) && ( dq > (median-rms * nsigma) ) ){
				truncated_dq += dq;
				npts += 1;
			}
		}

		dq_trunc_v.push_back( truncated_dq / npts );

		if(dq_trunc_v.back() != dq_trunc_v.back()){
			std::cout<<"ERROR::TruncMean.cxx || NAN "<<dq_trunc_v.back()<<std::endl;
			std::cout<<"truncated_dq: "<<truncated_dq<<std::endl;
			std::cout<<"npts: "<<npts<<std::endl;
			std::cout<<"median: "<<median<<std::endl;
			std::cout<<"rms: "<<rms<<std::endl;
		}

	}// for all values

	return;
}

double TruncMean::Mean(const std::vector<double>& v)
{

	double mean = 0.;
	for (auto const& n : v) mean += n;
	mean /= v.size();

	return mean;
}

double TruncMean::Median(const std::vector<double>& v)
{

	if (v.size() == 1) return v[0];

	std::vector<double> vcpy = v;

	std::sort(vcpy.begin(), vcpy.end());

	double median = vcpy[ vcpy.size() / 2 ];

	return median;
}

double TruncMean::RMS(const std::vector<double>& v)
{

	if(v.size()==1) return v.front();

	double avg = 0.;
	for (auto const& val : v) avg += val;
	avg /= v.size();
	double rms = 0.;
	for (auto const& val : v) rms += (val-avg)*(val-avg);
	rms = sqrt( rms / ( v.size() -  1 ) );


	if(rms!=rms){
		std::cout<<"ERROR || TruncMean::RMS || is returning nan."<<std::endl;
		std::cout<<"ERROR || TruncMean::RMS || "<<rms<<std::endl;
		std::cout<<"ERROR || TruncMean::RMS || Input is of size "<<v.size()<<std::endl;
		for(auto const &val : v){
			std::cout<<"ERROR || TruncMean::RMS || "<<val<<std::endl;
		}
		std::cout<<"ERROR || TruncMean::RMS || Most likely because your radius is too small!  "<<v.size()<<std::endl;
	}

	return rms;
}
#endif
/** @} */ // end of doxygen group 
