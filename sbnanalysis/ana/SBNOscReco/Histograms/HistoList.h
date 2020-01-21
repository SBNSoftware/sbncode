#ifndef _sbncode_histolisto_hh_
#define _sbncode_histolisto_hh_

#include <vector>

class TH1;
class TDirectory;

namespace ana {
 namespace SBNOsc {

class HistoList {
public:

  /**
  * Scale all histograms by a set value
  * \param scale The scaing value
  */
  void Scale(double scale);
  /**
  * Add another set of histograms to this one
  * \param other The set of histograms to add
  */
  void Add(const HistoList &other);
  
  /**
  ** Write this set of histograms to disk
  **/
  void Write();

  void StoreHisto(TH1 *histo);
  void Merge(const HistoList &merge);

  std::vector<TH1 *> fAllHistos;
  std::vector<TDirectory *> fLocations;
};

  } // end namespace ana
} // end namespace SBNOsc

#endif
