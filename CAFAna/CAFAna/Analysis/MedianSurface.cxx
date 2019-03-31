#include "CAFAna/Analysis/MedianSurface.h"

#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"

namespace ana
{
  // --------------------------------------------------------------------------
  MedianSurface::MedianSurface(const std::vector<Surface>& throws)
    : fThrows(throws)
  {
    assert(!throws.empty());

    for(const Surface& s: throws){
      assert(s.fHist->GetNbinsX() == throws[0].fHist->GetNbinsX());
      assert(s.fHist->GetNbinsY() == throws[0].fHist->GetNbinsY());
      // TODO check min and max at least
    }

    fLogX = throws[0].fLogX;
    fLogY = throws[0].fLogY;

    // Think about what to do with these
    fMinX = -1;
    fMinY = -1;
    // Though this is right
    fMinChi = 0;

    fHist = new TH2F(*throws[0].fHist);
    fHist->Reset();

    for(int ix = 0; ix < fHist->GetNbinsX()+2; ++ix){
      for(int iy = 0; iy < fHist->GetNbinsY()+2; ++iy){
        std::vector<float> chis;
        chis.reserve(throws.size());
        for(const Surface& s: throws){
          chis.push_back(s.fHist->GetBinContent(ix, iy));
        }
        std::sort(chis.begin(), chis.end());

        // Median
        fHist->SetBinContent(ix, iy, chis[throws.size()/2]);
      } // end for iy
    } // end for ix
  }

  // --------------------------------------------------------------------------
  void MedianSurface::DrawEnsemble(TH2* fc, Color_t color)
  {
    EnsureAxes();

    for(Surface& s: fThrows){
      std::vector<TGraph*> gs = s.GetGraphs(fc, -1);

      for(TGraph* g: gs){
        g->SetLineWidth(1);
        g->SetLineColor(color);
        g->Draw("l");
      }
    }

    gPad->Update();
  }
}
