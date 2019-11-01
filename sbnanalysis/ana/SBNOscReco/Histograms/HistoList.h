#ifndef _sbncode_histolisto_hh_
#define _sbncode_histolisto_hh_
class TH1;
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

  std::vector<TH1 *> fAllHistos;
};

  } // end namespace ana
} // end namespace SBNOsc

#endif
