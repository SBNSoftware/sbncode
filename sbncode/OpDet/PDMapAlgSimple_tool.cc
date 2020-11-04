///////////////////////////////////////////////////////////////////////
///
/// \file   PDMapAlg.h
///
/// \brief  This is the interface class for a tool to handle PD mapping
///         in SBN detectors, used for flash matching.
///
/// \author Laura Paulucci, Franciole Marinho, Iker de Icaza, and W. Ketchum
///
////////////////////////////////////////////////////////////////////////

#include "PDMapAlg.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include <string>

namespace opdet {

  //a default dummy one
  class PDMapAlgSimple : PDMapAlg
  {
  public:
    explicit PDMapAlgSimple(const fhicl::ParameterSet& pset)
    { fType = pset.get<std::string>("Type","pmt"); }
    
    ~PDMapAlgSimple(){}

    std::string pdType(size_t ch) const override { return fType; }

  private:
    std::string fType;
    
  };


  DEFINE_ART_CLASS_TOOL(PDMapAlgSimple)
} // namespace

