/**
 *  @file   IMeVPrtlDecay.h
 *
 *  @brief  This is an interface for an art Tool which decays "Prtl" inside a
 *  detector volume. It maps MeVPrtlFlux to MeVPrtlDecay. 
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IMeVPrtlDecay_h
#define IMeVPrtlDecay_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlDecay.h"

#include "IMeVPrtlStage.h"
#include "Constants.h"

// Algorithm includes

#include <utility>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {

/**
 *  @brief  IMeVPrtlDecay interface class definiton
 */
class IMeVPrtlDecay: virtual public IMeVPrtlStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IMeVPrtlDecay() noexcept = default;

    virtual bool Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) = 0;

protected:
    double TimeOfFlight(const MeVPrtlFlux &flux, TVector3 decay);
};

double IMeVPrtlDecay::TimeOfFlight(const MeVPrtlFlux &flux, TVector3 decay){

  // TODO: should the neutrino TOF be subtracted here to get the correct T0?
  // The current temporary fix in GENIE doesn't subtract first neutrino arrival time 
 
  //Work in detector coordinate
  //Calculate kaon ToF
  double kaonvel = flux.kpmom_beamcoord.Beta() * Constants::Instance().c_cm_per_ns;  
  TVector3 kaondecay = TVector3(flux.pos.X(), flux.pos.Y(), flux.pos.Z());  
  double kaonToF = kaondecay.Mag() / kaonvel;

  //Calculate distance of from flux to decay pos
  double hnlvel = flux.mom.Beta() * Constants::Instance().c_cm_per_ns; 
  double hnlToF = (decay - kaondecay).Mag() / hnlvel; // s -> ns

  std::cout << "=-----------------------------------------------------=" << std::endl;
  std::cout << "kaon ToF: " << kaonToF << std::endl;
  std::cout << "kaon vel: " << kaonvel << std::endl;
  std::cout << "=-----------------------------------------------------=" << std::endl;
  std::cout << "hnl ToF: " << hnlToF << std::endl;
  std::cout << "hnl vel: " << hnlvel << std::endl;
  std::cout << "=-----------------------------------------------------=" << std::endl;

  return kaonToF + hnlToF;
}//Time Of Flight



} // namespace ldm
} // namespace evgen
#endif

