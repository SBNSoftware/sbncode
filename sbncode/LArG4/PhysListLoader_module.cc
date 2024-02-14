////////////////////////////////////////////////////////////////////////
// Class:       PhysListLoader
// Plugin Type: producer (art v3_06_03)
// File:        PhysListLoader_module.cc
//
// Generated at Tue Aug  3 12:56:30 2021 by Marco Del Tutto using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

/**
 * @file PhysListLoader_module.cc
 *
 * @brief A dummy module that forces library loading to register a
 * physics list.
 *
 * @author Marco Del Tutto (mdeltutt@fnal.gov)
 *
 * This module is a dummy module that does not act in any way
 * on art::Events. This module, and in particular the library this
 * module is built in, is linked to the PhysicsLists library in
 * this same directory. Running this module causes the PhysicsLists
 * library to be loaded, and so the Physics Lists to be registered
 * to the Geant4 PhysListRegistry
 *
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "Geant4/G4PhysListRegistry.hh"


class PhysListLoader;


class PhysListLoader : public art::EDProducer {
public:
  explicit PhysListLoader(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhysListLoader(PhysListLoader const&) = delete;
  PhysListLoader(PhysListLoader&&) = delete;
  PhysListLoader& operator=(PhysListLoader const&) = delete;
  PhysListLoader& operator=(PhysListLoader&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


PhysListLoader::PhysListLoader(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.


  if (p.get<bool>("DumpList", false)) {
    std::cout << "[PhysListLoader] Dumping G4 physics list:" << std::endl;
    G4PhysListRegistry* g4plr = G4PhysListRegistry::Instance();
    g4plr->PrintAvailablePhysLists();
  }

}

void PhysListLoader::produce(art::Event& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(PhysListLoader)
