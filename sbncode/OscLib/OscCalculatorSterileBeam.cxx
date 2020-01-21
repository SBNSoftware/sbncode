

#include "OscCalculatorSterile.h"
#include "OscCalculatorSterileBeam.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <typeinfo>

namespace osc
{
  OscCalculatorSterileBeam::OscCalculatorSterileBeam()
    : OscCalculatorSterile(), 
    fKaonscale(0), fPionscale(0), fMuonscale(0)
  {
    //fPMNS_Sterile = new PMNS_Sterile(fNFlavors);
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterileBeam::~OscCalculatorSterileBeam()
  {
	  std::cout << " ** OscCalculatorSterileBeam Destructor is called ** \n";
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterileBeam::OscCalculatorSterileBeam(const OscCalculatorSterileBeam& calc)
    : OscCalculatorSterile(calc)
  {
    std::cout << "copy constructor\n";
    this->fKaonscale = calc.GetKaonScale();
    this->fMuonscale = calc.GetMuonScale();
    this->fPionscale = calc.GetPionScale();
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterileBeam::SetKaonScale(double scale)
  {
    fKaonscale=scale;
    //assert(fKaonscale!=0.0);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterileBeam::SetPionScale(double scale)
  {
    fPionscale=scale;
    //assert(fPionscale!=0.0);
  }

  //---------------------------------------------------------------------------
  void OscCalculatorSterileBeam::SetMuonScale(double scale)
  {
    fMuonscale=scale;
    //assert(fMuonscale!=0.0);
  }

  //---------------------------------------------------------------------------
  double OscCalculatorSterileBeam::GetKaonScale() const
  {
    return fKaonscale;
  }

  //---------------------------------------------------------------------------
  double OscCalculatorSterileBeam::GetPionScale() const
  {
    return fPionscale;
  }

  //---------------------------------------------------------------------------
  double OscCalculatorSterileBeam::GetMuonScale() const
  {
    return fMuonscale;
  }

  //---------------------------------------------------------------------------
  IOscCalculatorAdjustable* OscCalculatorSterileBeam::Copy() const
  {
    std::cout << " *** OscCalSterileBeam Copy() is called *** " << std::endl;
    return new OscCalculatorSterileBeam(*this);
  }

  //---------------------------------------------------------------------------
  TMD5* OscCalculatorSterileBeam::GetParamsHash() const
  {
    TMD5* ret = new TMD5;                                                            
    std::string txt = "SterileBeam";                                                 
    ret->Update((unsigned char*)txt.c_str(), txt.size());                            
    std::vector<double> buf;                                                         
    buf.push_back(fL);                                                               
    buf.push_back(fRho);                                                             
    buf.push_back(fKaonscale);                                                       
    buf.push_back(fPionscale);                                                       
    buf.push_back(fMuonscale);                                                       
    for(int i = 2; i <= fNFlavors; ++i) buf.push_back(GetDm(i));                     
    for(int j = 2; j <= fNFlavors; ++j) {                                            
      for(int i = 1; i < j; ++i) {                                                   
	buf.push_back(GetAngle(i, j));                                             
	if(i+1 != j) buf.push_back(GetDelta(i, j));                                
      }
    }
    ret->Update((unsigned char*)&buf[0], sizeof(double)*buf.size());                 
    ret->Final();                                                                    
    return ret;                                                                      
  }
  //---------------------------------------------------------------------------
  /* virtual void OscCalculatorSterileBeam::SetBeamMode(kBeamMode, double scaleK, double scaleP, double scaleM)
  {
    if(kBeamMode=="Kaon") {
      OscCalculatorSterileBeam::SetKaonScale(scaleK) ; 
      OscCalculatorSterileBeam::SetPionScale(1.0) ; 
      OscCalculatorSterileBeam::SetMuonScale(1.0) ;  }
    else if(kBeamMode=="Pion") {
      OscCalculatorSterileBeam::SetKaonScale(1.0) ; 
      OscCalculatorSterileBeam::SetPionScale(scaleP); 
      OscCalculatorSterileBeam::SetMuonScale(1.0) ;  }
    else if(kBeamMode=="Muon") {
      OscCalculatorSterileBeam::SetKaonScale(1.0) ; 
      OscCalculatorSterileBeam::SetPionScale(1.0) ; 
      OscCalculatorSterileBeam::SetMuonScale(scaleM); }
    else {
      OscCalculatorSterileBeam::SetKaonScale(1.0) ;  
      OscCalculatorSterileBeam::SetKaonScale(1.0) ;
      OscCalculatorSterileBeam::SetKaonScale(1.0) ;  }
  }
  
  //---------------------------------------------------------------------------
  virtual double OscCalculatorSterileBeam::P(int flavBefore, int flavAfter, double E)
  {
    double origProb = OscCalculatorSterile::P(int flavBefore, int flavAfter, double E);
    double scale = 0.0;
    double prob=0.0;
    if(kBeamMode=="Kaon" || kBeamMode=="Pion" || kBeamMode=="Muon") {
      scale = OscCalculatorSterileBeam::GetKaonScale() * 
	OscCalculatorSterileBeam::GetPionScale() * 
	OscCalculatorSterileBeam::GetMuonScale() ;
      prob=origProb*scale;    
    }
    return prob;
  }
  */
  //---------------------------------------------------------------------------
  const OscCalculatorSterileBeam* DowncastToSterileBeam(const IOscCalculator* calc)
  {
    //const IOscCalculator* base = new osc::OscCalculatorSterileBeam();
    const OscCalculatorSterileBeam* calc_sterile = dynamic_cast<const OscCalculatorSterileBeam*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalculatorSterileBeam." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

  //---------------------------------------------------------------------------
  OscCalculatorSterileBeam* DowncastToSterileBeam(IOscCalculator* calc)
  {
    //IOscCalculator* base = new osc::OscCalculatorSterileBeam();
    //std::cout << "calc's type is: " << typeid(calc).name() << std::endl;
    //if (typeid(calc) == typeid(osc::NoOscillations))
    //  {
    ///	std::cout << "calc is a osc::NoOscillations calculator!" << std::endl;
    //}
    NoOscillations* calc_noosc = dynamic_cast<NoOscillations*>(calc);
    if (calc_noosc) 
      {
	std::cout << "I was successfully cast to NoOscillations" << std::endl;
      }

    OscCalculatorSterileBeam* calc_sterile = dynamic_cast<OscCalculatorSterileBeam*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalculatorSterileBeam." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }
} // namespace
