// Check header file from documentation

#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include "Geant4/globals.hh"
#include "Geant4/G4ios.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4ProcessVector.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ParticleTable.hh"

#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialTable.hh"

#include "Geant4/G4DecayPhysics.hh"
#include "Geant4/G4EmStandardPhysics.hh"
#include "Geant4/G4EmExtraPhysics.hh"
#include "Geant4/G4IonPhysics.hh"
#include "Geant4/G4StoppingPhysics.hh"
#include "Geant4/G4HadronElasticPhysicsHP.hh"
// #include "Geant4/G4NeutronTrackingCut.hh"
#include "Geant4/G4HadronPhysicsQGSP_BERT_HP.hh"


// Register this new physics list
#include "Geant4/G4PhysListStamper.hh" // defines macro for factory registration
#include "SBN_QGSP_BERT_NNC.hh"
G4_DECLARE_PHYSLIST_FACTORY(SBN_QGSP_BERT_NNC);

SBN_QGSP_BERT_NNC::SBN_QGSP_BERT_NNC(G4int ver)
{
  if(ver > 0) {
    G4cout << "<<< Geant4 Physics List simulation engine: SBN_QGSP_BERT_NNC"<<G4endl;
    G4cout <<G4endl;
  }

  defaultCutValue = 0.7*CLHEP::mm;
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );

  // Hadron Physics
  RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  // Neutron tracking cut
  // RegisterPhysics( new G4NeutronTrackingCut(ver));

}
