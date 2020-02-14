{
  // event_numbers.C has to be executed first!
  // Replace SBND with the appropriate experiment
  TFile* f = new TFile("output/output_finalstate_SBND.root");
  
  TH1D* SBND_All_NC_all = (TH1D*) f->Get("SBND_All_NC_all");
  TH1D* SBND_0pi_CC_all = (TH1D*) f->Get("SBND_0pi_CC_all");
  TH1D* SBND_1piCh_CC_all = (TH1D*) f->Get("SBND_1piCh_CC_all");
  TH1D* SBND_2piCh_CC_all = (TH1D*) f->Get("SBND_2piCh_CC_all");
  TH1D* SBND_1each_CC_all = (TH1D*) f->Get("SBND_1each_CC_all");
  TH1D* SBND_1pi0_CC_all = (TH1D*) f->Get("SBND_1pi0_CC_all");
  TH1D* SBND_2pi0_CC_all = (TH1D*) f->Get("SBND_2pi0_CC_all");
  TH1D* SBND_else_CC_all = (TH1D*) f->Get("SBND_else_CC_all");
  TH1D* SBND_All_CC_all = (TH1D*) f->Get("SBND_All_CC_all");

  SBND_All_NC_all->SetFillColorAlpha(SBND_All_NC_all->GetLineColor(), 0.6);
  SBND_0pi_CC_all->SetFillColorAlpha(SBND_0pi_CC_all->GetLineColor(), 0.6);
  SBND_1piCh_CC_all->SetFillColorAlpha(SBND_1piCh_CC_all->GetLineColor(), 0.6);
  SBND_2piCh_CC_all->SetFillColorAlpha(SBND_2piCh_CC_all->GetLineColor(), 0.6);
  SBND_1each_CC_all->SetFillColorAlpha(SBND_1each_CC_all->GetLineColor(), 0.6);
  SBND_1pi0_CC_all->SetFillColorAlpha(SBND_1pi0_CC_all->GetLineColor(), 0.6);
  SBND_2pi0_CC_all->SetFillColorAlpha(SBND_2pi0_CC_all->GetLineColor(), 0.6);
  SBND_else_CC_all->SetFillColorAlpha(SBND_else_CC_all->GetLineColor(), 0.6);

  //Change POT if needed
  THStack* h = new THStack("h",";Reconstructed energy (GeV);Events / (6.6#times10^{20} POT)");

  h->Add(SBND_All_NC_all);
  h->Add(SBND_0pi_CC_all);
  h->Add(SBND_1piCh_CC_all);
  h->Add(SBND_2piCh_CC_all);
  h->Add(SBND_1pi0_CC_all);
  h->Add(SBND_2pi0_CC_all);
  h->Add(SBND_1each_CC_all);
  h->Add(SBND_else_CC_all);

  TLegend ll(0.65,0.5,0.84,0.9);

  double n = SBND_All_CC_all->Integral() + SBND_All_NC_all->Integral();
  double nNC = SBND_All_NC_all->Integral();
  double n0pi = SBND_0pi_CC_all->Integral();
  double n1piCh = SBND_1piCh_CC_all->Integral();
  double n2piCh = SBND_2piCh_CC_all->Integral();
  double n1pi0 = SBND_1pi0_CC_all->Integral();
  double n2pi0 = SBND_2pi0_CC_all->Integral();
  double n1each = SBND_1each_CC_all->Integral();
  double nelse = SBND_else_CC_all->Integral();

  ll.SetFillColor(10);
  ll.AddEntry(SBND_All_NC_all, Form("NC (%3.1f%%)",nNC*100./n), "f");
  ll.AddEntry(SBND_0pi_CC_all, Form("CC 0#pi (%3.1f%%)",n0pi*100./n), "f");
  ll.AddEntry(SBND_1piCh_CC_all, Form("CC 1#pi^{#pm} (%3.1f%%)",n1piCh*100./n), "f");
  ll.AddEntry(SBND_2piCh_CC_all, Form("CC 2#pi^{#pm} (%3.1f%%)",n2piCh*100./n), "f");
  ll.AddEntry(SBND_1pi0_CC_all, Form("CC 1#pi^{0} (%3.1f%%)",n1pi0*100./n), "f");
  ll.AddEntry(SBND_2pi0_CC_all, Form("CC 2#pi^{0} (%3.1f%%)",n2pi0*100./n), "f");
  ll.AddEntry(SBND_1each_CC_all, Form("CC 1#pi^{0}1#pi^{#pm} (%3.1f%%)",n1each*100./n), "f");
  ll.AddEntry(SBND_else_CC_all, Form("CC else (%3.1f%%)",nelse*100./n), "f");

  TLatex* prelim = new TLatex(.9, .96, "SBND Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);

  TCanvas* c1 = new TCanvas("c1");
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();
  c1->SaveAs("output/Stack_finalstate_numu.root");
  c1->SaveAs("output/Stack_finalstate_numu.pdf");
  
  TCanvas* c2 = new TCanvas("c2","c2", 700,800);
  c2->Divide(1,2);

  c2->cd(1);
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();

  c2->cd(2);
  h->Draw("nostack hist");
  prelim->Draw();
  c2->SaveAs("output/StackNostack_finalstate_numu.pdf");
  c2->SaveAs("output/StackNostack_finalstate_numu.pdf");
}
