//make_FD_reco_systs.C
// macro to make the histograms for the reconstruction systematics based upon differences in 
// kinematic variables between generators

// Makes most of the histograms
void make_hists(TTree* inputTree, const char* treeName) {
  // 1D hists
  TH1D *hx  = new TH1D(Form("hx_%s", treeName), Form("X distribution for %s; X; Events", treeName), 40, 0, 2.5);
  TH1D *hy  = new TH1D(Form("hy_%s", treeName), Form("Y distribution for %s; Y; Events", treeName), 20, 0, 1);
  TH1D *hw  = new TH1D(Form("hw_%s", treeName), Form("W distribution for %s; W; Events", treeName), 40, 0, 4);
  TH1D *hq2 = new TH1D(Form("hq2_%s", treeName), Form("Q^{2} distribution for %s; Q^{2}; Events", treeName), 20, 0, 8);
  hx->SetOption("hist");
  hy->SetOption("hist");
  hw->SetOption("hist");
  hq2->SetOption("hist");
  hx->SetStats(kFALSE);
  hy->SetStats(kFALSE);
  hw->SetStats(kFALSE);
  hq2->SetStats(kFALSE);
  // 2D hists
  TH2D *hEvX = new TH2D(Form("hEvX_%s", treeName), Form("X vs E_{#nu} for %s, normalized to E_{#nu}; E_{#nu}; X", treeName), 20, 0, 10, 20, 0, 2.5);
  TH2D *hEvY = new TH2D(Form("hEvY_%s", treeName), Form("Y vs E_{#nu} for %s, normalized to E_{#nu}; E_{#nu}; Y", treeName), 20, 0, 10, 20, 0, 1);
  TH2D *hEvW = new TH2D(Form("hEvW_%s", treeName), Form("W vs E_{#nu} for %s, normalized to E_{#nu}; E_{#nu}; W", treeName), 20, 0, 10, 20, 0, 4);
  TH2D *hEvQ2 = new TH2D(Form("hEvQ2_%s", treeName), Form("Q^{2} vs E_{#nu} for %s, normalized to E_{#nu}; E_{#nu}; Q^{2}", treeName), 20, 0, 10, 20, 0, 7);
  hEvX->SetOption("colz");
  hEvY->SetOption("colz");
  hEvW->SetOption("colz");
  hEvQ2->SetOption("colz");
  hEvX->SetStats(kFALSE);
  hEvY->SetStats(kFALSE);
  hEvQ2->SetStats(kFALSE);
  hEvW->SetStats(kFALSE);

  double Enu;
  int nfsp;
  int pdg[100];
  double px[100];
  double py[100];
  double pz[100];
  double E[100];
  inputTree->SetBranchAddress("Enu", &Enu);
  inputTree->SetBranchAddress("nfsp", &nfsp);
  inputTree->SetBranchAddress("pdg", pdg);
  inputTree->SetBranchAddress("px", px);
  inputTree->SetBranchAddress("py", py);
  inputTree->SetBranchAddress("pz", pz);
  inputTree->SetBranchAddress("E", E);

  const TLorentzVector p_p(0, 0, 0, .9315); // proton at rest

  std::cout<<"Calculating kinematics for "<<treeName<<std::endl;
  for (int t=0; t<inputTree->GetEntries(); t++) {
    inputTree->GetEntry(t);
    TLorentzVector lep(0., 0., 0., 0.); 
    TLorentzVector nu;
    nu.SetE(Enu);
    nu.SetPz(Enu);
    double Sumpx = 0.;
    double Sumpy = 0.;
    double Sumpz = 0.;
    double SumE  = 0.;
    bool isCC = false;
    for (int i=0; i<nfsp; i++) {
      if (abs(pdg[i]) == 11 || abs(pdg[i]) == 13 || abs(pdg[i]) == 15) {
	lep.SetPxPyPzE(px[i], py[i], pz[i], E[i]);
	isCC = true;
      }
      else { }
    }
    // Calculate kinematics
    if (isCC) {
      const TLorentzVector q = nu - lep;
      const double Qsq = -(q.Mag2());
      const double w = (q + p_p).Mag();
      const double x = Qsq / (2*p_p.Dot(q));
      const double y = 1 - lep.E()/nu.E();

      hq2->Fill(Qsq);
      hw->Fill(w);
      hx->Fill(x);
      hy->Fill(y);
      hEvQ2->Fill(nu.E(), Qsq);
      hEvW->Fill(nu.E(), w);
      hEvX->Fill(nu.E(), x);
      hEvY->Fill(nu.E(), y);
    }
  }
  // Normalize
  double integral = 1. / hq2->Integral();
  hq2->Scale(integral);
  hw->Scale(integral);
  hx->Scale(integral);
  hy->Scale(integral);
  for (int col=0; col <= hEvQ2->GetNbinsX(); col++) {
    double EvInt = hEvQ2->Integral(col, col, 0, 20);
    for (int row=0; row <= hEvQ2->GetNbinsY(); row++) {
      hEvQ2->SetBinContent(col, row, hEvQ2->GetBinContent(col, row) / EvInt);
      hEvW->SetBinContent(col, row, hEvW->GetBinContent(col, row) / EvInt);
      hEvX->SetBinContent(col, row, hEvX->GetBinContent(col, row) / EvInt);
      hEvY->SetBinContent(col, row, hEvY->GetBinContent(col, row) / EvInt);
    }
  }

  hq2->Write();
  hw->Write();
  hx->Write();
  hy->Write();
  hEvQ2->Write();
  hEvW->Write();
  hEvX->Write();
  hEvY->Write();
}
// Takes 2 2D spectra and makes the corresponding ratio plot
void make_ratio_plot(TFile* f, const char* nomname, const char* denomname, const char* var, const char* nom, const char* denom) {
  TH2 *hnom   = (TH2*)f->Get(nomname);
  TH2 *hdenom = (TH2*)f->Get(denomname);
  int xBins = hnom->GetNbinsX();
  int yBins = hnom->GetNbinsY();
  double xMax  = hnom->GetXaxis()->GetXmax();
  double yMax  = hnom->GetYaxis()->GetXmax();

  TH2D *hratio = new TH2D(Form("h%sratio_%s_%s", var, nom, denom), Form("%s/%s ratio for %s", nom, denom, var), xBins, 0, xMax, yBins, 0, yMax);
  hratio->SetOption("colz");
  hratio->SetStats(kFALSE);
  hratio->Divide(hnom, hdenom);
  hratio->Write();
}
// Makes all the relevant systematics plots
// Current location /pnfs/dune/persistent/users/marshalc/modelComp/
void make_FD_reco_systs(const char* geniefile, const char* neutfile) {
  TFile *gfile = new TFile(geniefile, "read");
  TFile *nfile = new TFile(neutfile, "read");
  TTree *neutTreefhc = (TTree*)nfile->Get("NEUT_FHC");
  TTree *neutTreerhc = (TTree*)nfile->Get("NEUT_RHC");
  TTree *genieTreefhc = (TTree*)gfile->Get("GENIE_FHC");
  TTree *genieTreerhc = (TTree*)gfile->Get("GENIE_RHC");

  TFile *fout = new TFile("modelComp.root", "recreate");

  make_hists(genieTreefhc, "geniefhc");
  make_hists(genieTreerhc, "genierhc");
  make_hists(neutTreefhc, "neutfhc");
  make_hists(neutTreerhc, "neutrhc");

  std::cout<<"Making ratio plots"<<std::endl;
  make_ratio_plot(fout, "hEvQ2_neutfhc", "hEvQ2_geniefhc", "Q2", "neutfhc", "geniefhc");
  make_ratio_plot(fout, "hEvQ2_neutrhc", "hEvQ2_genierhc", "Q2", "neutrhc", "genierhc");
  make_ratio_plot(fout, "hEvW_neutfhc", "hEvW_geniefhc", "W", "neutfhc", "geniefhc");
  make_ratio_plot(fout, "hEvW_neutrhc", "hEvW_genierhc", "W", "neutrhc", "genierhc");
  make_ratio_plot(fout, "hEvX_neutfhc", "hEvX_geniefhc", "X", "neutfhc", "geniefhc");
  make_ratio_plot(fout, "hEvX_neutrhc", "hEvX_genierhc", "X", "neutrhc", "genierhc");
  make_ratio_plot(fout, "hEvY_neutfhc", "hEvY_geniefhc", "Y", "neutfhc", "geniefhc");
  make_ratio_plot(fout, "hEvY_neutrhc", "hEvY_genierhc", "Y", "neutrhc", "genierhc");

  fout->Write();
}






