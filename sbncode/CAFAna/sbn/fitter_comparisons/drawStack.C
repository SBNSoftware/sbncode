{
  // event_numbers.C has to be executed first!
  // Replace SBND with the appropriate experiment

  TFile* f = new TFile("output/output_SBND.root");
  
  TH1D* SBND_All_NC_all = (TH1D*) f->Get("SBND_All_NC_all");
  TH1D* SBND_QE_CC_all = (TH1D*) f->Get("SBND_QE_CC_all");
  TH1D* SBND_Res_CC_all = (TH1D*) f->Get("SBND_Res_CC_all");
  TH1D* SBND_DIS_CC_all = (TH1D*) f->Get("SBND_DIS_CC_all");
  TH1D* SBND_Coh_CC_all = (TH1D*) f->Get("SBND_Coh_CC_all");
  TH1D* SBND_MEC_CC_all = (TH1D*) f->Get("SBND_MEC_CC_all");
  TH1D* SBND_All_CC_all = (TH1D*) f->Get("SBND_All_CC_all");

  SBND_All_NC_all->SetFillColorAlpha(SBND_All_NC_all->GetLineColor(), 0.6);
  SBND_QE_CC_all->SetFillColorAlpha(SBND_QE_CC_all->GetLineColor(), 0.6);
  SBND_Res_CC_all->SetFillColorAlpha(SBND_Res_CC_all->GetLineColor(), 0.6);
  SBND_DIS_CC_all->SetFillColorAlpha(SBND_DIS_CC_all->GetLineColor(), 0.6);
  SBND_Coh_CC_all->SetFillColorAlpha(SBND_Coh_CC_all->GetLineColor(), 0.6);
  SBND_MEC_CC_all->SetFillColorAlpha(SBND_MEC_CC_all->GetLineColor(), 0.6);

  //Change POT if needed
  THStack* h = new THStack("h",";Reconstructed energy (GeV);Events / (6.6#times10^{20} POT)");

  h->Add(SBND_All_NC_all);
  h->Add(SBND_QE_CC_all);
  h->Add(SBND_Res_CC_all);
  h->Add(SBND_MEC_CC_all);
  h->Add(SBND_DIS_CC_all);
  h->Add(SBND_Coh_CC_all);

  TLegend ll(0.65,0.5,0.84,0.9);

  double n = SBND_All_CC_all->Integral() + SBND_All_NC_all->Integral();
  double nNC = SBND_All_NC_all->Integral();
  double nQE = SBND_QE_CC_all->Integral();
  double nRes = SBND_Res_CC_all->Integral();
  double nDIS = SBND_DIS_CC_all->Integral();
  double nMEC = SBND_MEC_CC_all->Integral();
  double nCoh = SBND_Coh_CC_all->Integral();

  ll.SetFillColor(10);
  ll.AddEntry(SBND_All_NC_all, Form("NC (%3.1f%%)",nNC*100./n), "f");
  ll.AddEntry(SBND_QE_CC_all, Form("CC QE (%3.1f%%)",nQE*100./n), "f");
  ll.AddEntry(SBND_Res_CC_all, Form("CC Res (%3.1f%%)",nRes*100./n), "f");
  ll.AddEntry(SBND_MEC_CC_all, Form("CC MEC (%3.1f%%)",nMEC*100./n), "f");
  ll.AddEntry(SBND_DIS_CC_all, Form("CC DIS (%3.1f%%)",nDIS*100./n), "f");
  ll.AddEntry(SBND_Coh_CC_all, Form("CC Coh (%3.1f%%)",nCoh*100./n), "f");

  TLatex* prelim = new TLatex(.9, .96, "SBND Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);

  TCanvas* c1 = new TCanvas("c1");
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();
  c1->SaveAs("output/Stack_numu.root");
  c1->SaveAs("output/Stack_numu.pdf");
  
  TCanvas* c2 = new TCanvas("c2","c2", 700,800);
  c2->Divide(1,2);

  c2->cd(1);
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();

  c2->cd(2);
  h->Draw("nostack hist");
  prelim->Draw();
  c2->SaveAs("output/StackNostack_numu.pdf");
  c2->SaveAs("output/StackNostack_numu.pdf");
}
