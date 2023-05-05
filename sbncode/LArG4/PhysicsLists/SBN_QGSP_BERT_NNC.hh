/**
 * @file SBN_QGSP_BERT_NNC
 *
 * @brief A Geant4 physics list. Same as QGSP_BERT,
 * but without neutron cut.
 *
 * @author Marco Del Tutto (mdeltutt@fnal.gov)
 *
 * A Geant4 physics list based on:
 * geant4/source/physics_lists/lists/src/QGSP_BERT.*
 * and removing the G4NeutronTrackingCut.
 *
 */

#ifndef SBN_QGSP_BERT_NNC_h
#define SBN_QGSP_BERT_NNC_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4VModularPhysicsList.hh"

class SBN_QGSP_BERT_NNC: public G4VModularPhysicsList
{
public:
  SBN_QGSP_BERT_NNC(G4int ver = 1);
  virtual ~SBN_QGSP_BERT_NNC()=default;

  SBN_QGSP_BERT_NNC(const SBN_QGSP_BERT_NNC &) = delete;
  SBN_QGSP_BERT_NNC & operator=(const SBN_QGSP_BERT_NNC &)=delete;

};

#endif
