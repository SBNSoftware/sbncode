////////////////////////////////////////////////////////////////////////
// \author  Bruno Zamorano
// \date    February 2018
////////////////////////////////////////////////////////////////////////
#ifndef SRNEUTRINO_H
#define SRNEUTRINO_H

namespace caf
{
  /// The SRNeutrino is a representation of neutrino interaction information
  class SRNeutrino
    {
    public:
      SRNeutrino();
      ~SRNeutrino() {  };

      bool         iscc;          ///< Is CC interaction
      int          pdg;           ///< PDG code
      int          genie_intcode; ///< GENIE interaction mode
      double       energy;        ///< True energy [GeV]
      double       inelasticityY; ///< True inelasticity
    };

} // end namespace

#endif // SRNEUTRINO_H
//////////////////////////////////////////////////////////////////////////////
