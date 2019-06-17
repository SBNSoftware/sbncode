#pragma once

#include <map>
#include <vector>

#include "Rtypes.h"

class TLegend;
class TGraph;
class TGraphAsymmErrors;
class TH1;
class TH2;
class THStack;

namespace osc{class IOscCalculator;}

namespace ana
{
  class IPrediction;
  class ISyst;
  class Spectrum;
  class SystShifts;

  /// Overlay MC spectrum with data spectrum
  ///
  /// \return The first histogram drawn so you can alter axis labels etc
  TH1* DataMCComparison(const Spectrum& data, const Spectrum& mc);

  /// Overlay MC spectrum with data spectrum, area normalized
  ///
  /// \return The first histogram drawn so you can alter axis labels etc
  TH1* DataMCComparisonAreaNormalized(const Spectrum& data, const Spectrum& mc);

  /// \brief Plot MC broken down into flavour components, overlayed with data
  ///
  /// \return The first histogram drawn so you can alter axis labels etc
  TH1* DataMCComparisonComponents(const Spectrum& data,
                                  const IPrediction* mc,
                                  osc::IOscCalculator* calc);

  TH1* GetMCSystTotal(const IPrediction* mc,
                      osc::IOscCalculator* calc,
                      const SystShifts& shift,
		      std::string hist_name,
                      double pot,
                      bool force1D=false);

  TH1* GetMCTotal(const IPrediction* mc,
                  osc::IOscCalculator* calc,
                  std::string hist_name,
                  double pot,
                  bool force1D=false);


  /// A vector of histograms for the MC components. Easier to manipulate elsewhere
  /// Not ideal as returned pointers probably won't get deleted... but very useful...
  std::vector<TH1*> GetMCComponents(const IPrediction* mc,
				    osc::IOscCalculator* calc,
				    std::string hist_name,
				    double pot,
				    bool force1D = false);
  
  std::vector<TH1*> GetMCSystComponents(const IPrediction* mc,
					osc::IOscCalculator* calc,
					const SystShifts& shift,
					std::string hist_name,
					double pot,
					bool force1D=false);


  std::vector<TH1*> GetMCTotalForSystShifts(const IPrediction* mc,
					    osc::IOscCalculator* calc,
					    const ISyst* syst,
					    std::string hist_base_name,
					    double pot,
					    bool force1D = false);

  /// Plot data/MC ratio for the given spectrum. Normalize MC to Data by area
  void DataMCAreaNormalizedRatio(const Spectrum& data, const Spectrum& mc,
                                 double miny = 0, double maxy = 3);

  /// Plot data/MC ratio for the given spectrum. Normalize MC to Data by area
  void DataMCAreaNormalizedRatio(const Spectrum& data,
                                 const IPrediction* mc,
                                 osc::IOscCalculator* calc,
                                 double miny = 0, double maxy = 3);

  /// Plot data/MC ratio for the given spectrum. Normalize MC to Data by POT
  void DataMCRatio(const Spectrum& data, const Spectrum& mc,
                   double miny = 0, double maxy = 3);

  /// Plot data/MC ratio for the given spectrum. Normalize MC to Data by POT
  void DataMCRatio(const Spectrum& data,
		   const IPrediction* mc,
		   osc::IOscCalculator* calc,
                   double miny = 0, double maxy = 3);

  /// Plot data/expected, compared with fit/expected
  void RatioPlot(const Spectrum& data,
                 const Spectrum& expected,
                 const Spectrum& fit,
                 double miny = 0, double maxy = 1.2);

  /// \brief Plot prediction with +/-1sigma error band
  ///
  /// When multiple systematics are used, the errors are the quadrature sum
  ///
  /// \param pred   The prediction. Must be able to generate shifted spectra
  /// \param systs  List of systematics to include in error band
  /// \param calc   Oscillation parameters
  /// \param pot    POT to evaluate prediction at
  /// \param col    Color of the prediction, default kRed
  /// \param errCol Color of the shading, default col-7 (kRed-7 is light red)
  /// \param headroom   Fraction of maximum bin for headroom, default 30%  
  void PlotWithSystErrorBand(IPrediction* pred,
                             const std::vector<const ISyst*>& systs,
                             osc::IOscCalculator* calc,
                             double pot,
                             int col = -1, int errCol = -1, 
                             float headroom = 1.3,
			     bool newaxis = true);

  /// \brief Plot prediction with error band
  ///
  /// When multiple systematics are used, the errors are the quadrature sum
  ///
  /// \param nominal    Nominal spectrum 
  /// \param upShifts   Vector of spectra which have + shifts  
  /// \param downShifts Vector of spectra which have - shifts, same order as +
  /// \param pot        POT to scale spectra to 
  /// \param col        Color of the prediction, default kRed
  /// \param errCol     Color of the shading, default col-7(kRed-7 is light red)
  /// \param headroom   Fraction of maximum bin for headroom, default 30%
  void PlotWithSystErrorBand(const Spectrum& nominal, 
                             const std::vector<Spectrum>& upShifts,
                             const std::vector<Spectrum>& downShifts,
                             double pot,
                             int col = -1, int errCol = -1,
                             float headroom = 1.3, bool newaxis=true);

  /// Can call like ToTHStack({{h1, kRed}, {h2, kBlue}}, pot)
  THStack* ToTHStack(const std::vector<std::pair<Spectrum, Color_t>>& s,
                     double pot);


  /// \brief Create a legend, maximizing distance from all histograms
  ///
  /// \param dx Width of the legend, fraction of pad width
  /// \param dy Height of the legend, fraction of pad height
  /// \param yPin, height in NDC (frac) to pin center of legend, i.e only move x
  TLegend* AutoPlaceLegend(double dx, double dy, double yPin = -1);

  /// \brief Make a simple plot of relative size of different errors
  ///
  /// \param systs     Map from error label to size in percent
  /// \param statErr   Percentage statistical error. Optional.
  /// \param bkgdOrSig Label axis for background or signal. Default (false) = bkgd
  void CountingExperimentErrorBarChart(const std::map<std::string, double>& systs, double statErr = 0, bool bkgdOrSig = false, bool shortchart=false);

  /// Calculate statistical errors appropriate for small Poisson numbers
  TGraphAsymmErrors* GraphWithPoissonErrors(const TH1* h, bool noErrorsXaxis = false, bool drawEmptyBins = true);

  /// Gives a TGraph with the area between two histograms. Do Draw("f") to draw 
  /// this area. By default it has a lighter version of the colour of hmin
  TGraph* ShadeBetweenHistograms(TH1* hmin, TH1* hmax);

  /// Calculate profile with error bars corresponding to specified quantiles of a 2D distribution (by default, 68% coverage)
  TGraphAsymmErrors * ProfileQuantile(const TH2 * hist, const std::string & axis_name, const std::string & graph_name="", const std::pair<double, double> & quantile_divisions={0.159, 0.841});
}
