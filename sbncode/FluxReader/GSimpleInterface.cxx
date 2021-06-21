#include "GSimpleInterface.h"
#include "TFile.h"
#include "TTree.h"

namespace fluxr {
  GSimpleInterface::GSimpleInterface()
  {
  }

  GSimpleInterface::~GSimpleInterface()
  {
  }

  void GSimpleInterface::SetRootFile(TFile* fluxInputFile)
  {
    fFluxTree=dynamic_cast<TTree*>(fluxInputFile->Get("flux"));
    fMetaTree=dynamic_cast<TTree*>(fluxInputFile->Get("meta"));
    fGSimpleEntry = new genie::flux::GSimpleNtpEntry;
    fGSimpleNuMI  = new genie::flux::GSimpleNtpNuMI;
    fGSimpleAux   = new genie::flux::GSimpleNtpAux;
    fGSimpleMeta  = new genie::flux::GSimpleNtpMeta;
    fFluxTree->SetBranchAddress("entry",&fGSimpleEntry);
    fFluxTree->SetBranchAddress("numi" ,&fGSimpleNuMI);
    fFluxTree->SetBranchAddress("aux"  ,&fGSimpleAux);
    fMetaTree->SetBranchAddress("meta" ,&fGSimpleMeta);

    fNEntries=fFluxTree->GetEntries();
    fFluxTree->GetEntry(0);
    fRun=fGSimpleNuMI->run;
    fMetaTree->GetEntry(0);
    fPOT=fGSimpleMeta->protons;

    Double_t* window_base = fGSimpleMeta->windowBase;
    Double_t* window_dir_1 = fGSimpleMeta->windowDir1;
    Double_t* window_dir_2 = fGSimpleMeta->windowDir2;
    std::cout << "GSimple Window:" << std::endl;
    std::cout << "\tBase: " << window_base[0] << ", " << window_base[1] << ", " << window_base[2] << std::endl;
    std::cout << "\tDir 1: " << window_dir_1[0] << ", " << window_dir_1[1] << ", " << window_dir_1[2] << std::endl;
    std::cout << "\tDir 2: " << window_dir_2[0] << ", " << window_dir_2[1] << ", " << window_dir_2[2] << std::endl;
  }

  bool GSimpleInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux)
  {
    if (!fFluxTree->GetEntry(ientry))
      return false;

    fNuPos=TLorentzVector(fGSimpleEntry->vtxx*100,fGSimpleEntry->vtxy*100,fGSimpleEntry->vtxz*100,0);
    fNuMom=TLorentzVector(fGSimpleEntry->px  ,fGSimpleEntry->py  ,fGSimpleEntry->pz  ,fGSimpleEntry->E);

    //stealing code from GenieHelper (nutools/EventGeneratorBase/GENIE/GENIEHelper.cxx)
    flux.fntype  = fGSimpleEntry->pdg;
    flux.fnimpwt = fGSimpleEntry->wgt;
    flux.fdk2gen = fGSimpleEntry->dist;
    flux.fnenergyn = flux.fnenergyf = fGSimpleEntry->E;

    if ( fGSimpleNuMI ) {
      flux.frun      = fGSimpleNuMI->run;
      flux.fevtno    = fGSimpleNuMI->evtno;
      flux.ftpx      = fGSimpleNuMI->tpx;
      flux.ftpy      = fGSimpleNuMI->tpy;
      flux.ftpz      = fGSimpleNuMI->tpz;
      flux.ftptype   = fGSimpleNuMI->tptype;   // converted to PDG
      flux.fvx       = fGSimpleNuMI->vx;
      flux.fvy       = fGSimpleNuMI->vy;
      flux.fvz       = fGSimpleNuMI->vz;

      flux.fndecay   = fGSimpleNuMI->ndecay;
      flux.fppmedium = fGSimpleNuMI->ppmedium;

      flux.fpdpx     = fGSimpleNuMI->pdpx;
      flux.fpdpy     = fGSimpleNuMI->pdpy;
      flux.fpdpz     = fGSimpleNuMI->pdpz;

      double apppz = fGSimpleNuMI->pppz;
      if ( TMath::Abs(fGSimpleNuMI->pppz) < 1.0e-30 ) apppz = 1.0e-30;
      flux.fppdxdz   = fGSimpleNuMI->pppx / apppz;
      flux.fppdydz   = fGSimpleNuMI->pppy / apppz;
      flux.fpppz     = fGSimpleNuMI->pppz;

      flux.fptype    = fGSimpleNuMI->ptype;
    }

    // anything useful stuffed into vdbl or vint?
    // need to check the metadata  auxintname, auxdblname
    if ( fGSimpleAux && fGSimpleMeta ) {
      // references just for reducing complexity
      const std::vector<std::string>& auxdblname = fGSimpleMeta->auxdblname;
      const std::vector<std::string>& auxintname = fGSimpleMeta->auxintname;
      const std::vector<int>&    auxint = fGSimpleAux->auxint;
      const std::vector<double>& auxdbl = fGSimpleAux->auxdbl;

      for (size_t id=0; id<auxdblname.size(); ++id) {
	if ("muparpx"   == auxdblname[id]) flux.fmuparpx  = auxdbl[id];
	if ("muparpy"   == auxdblname[id]) flux.fmuparpy  = auxdbl[id];
	if ("muparpz"   == auxdblname[id]) flux.fmuparpz  = auxdbl[id];
	if ("mupare"    == auxdblname[id]) flux.fmupare   = auxdbl[id];
	if ("necm"      == auxdblname[id]) flux.fnecm     = auxdbl[id];
	if ("nimpwt"    == auxdblname[id]) flux.fnimpwt   = auxdbl[id];
	if ("fgXYWgt"   == auxdblname[id]) {
	  flux.fnwtnear = flux.fnwtfar = auxdbl[id];
	}
      }
      for (size_t ii=0; ii<auxintname.size(); ++ii) {
	if ("tgen"      == auxintname[ii]) flux.ftgen     = auxint[ii];
	if ("tgptype"   == auxintname[ii]) flux.ftgptype  = auxint[ii];
      }
    }
    return true;
  }
}
