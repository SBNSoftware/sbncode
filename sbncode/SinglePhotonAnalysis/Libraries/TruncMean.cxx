#include "TruncMean.h"

namespace single_photon
{

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

    for (size_t n=0; n < dq_v.size(); n++) {   // rr_v and dq_v have the same size

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

        dq_local_v.push_back( dq_v[i] );  //save for ticks that are close enough

      }// for all ticks we want to scan

      if (dq_local_v.size() == 0 ) {
        dq_trunc_v.push_back( dq_v.at(n) ); // if no neighbours, push back dq of itself
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

      dq_trunc_v.push_back( truncated_dq / npts ); // push back averaged dq for these sitting within nsigma

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

    std::sort(vcpy.begin(), vcpy.end());  //sort to ascending order

    double median = vcpy[ vcpy.size() / 2 ]; // what if v has even # of elements? there is a choice

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

}
