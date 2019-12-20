{
  // event_numbers.C has to be executed first!
  // Replace SBND with the appropriate experiment

  TFile* f = new TFile("output/output_nue_SBND.root");
  
  TH1D* SBND_NC = (TH1D*) f->Get("SBND_NC");
  TH1D* SBND_numu = (TH1D*) f->Get("SBND_numu");
  TH1D* SBND_int = (TH1D*) f->Get("SBND_int");
  TH1D* SBND_osc = (TH1D*) f->Get("SBND_osc");

  SBND_NC->SetFillColorAlpha(SBND_NC->GetLineColor(), 0.6);
  SBND_numu->SetFillColorAlpha(SBND_numu->GetLineColor(), 0.6);
  SBND_int->SetFillColorAlpha(SBND_int->GetLineColor(), 0.6);
  SBND_osc->SetFillColorAlpha(SBND_osc->GetLineColor(), 0.6);

  //Change POT if needed
  THStack* h = new THStack("h",";Reconstructed energy (GeV);Events / (6.6#times10^{20} POT)");

  h->Add(SBND_int);
  h->Add(SBND_NC);
  h->Add(SBND_numu);
  h->Add(SBND_osc);

  TLegend ll(0.65,0.5,0.84,0.9);

  double n = SBND_NC->Integral() + SBND_numu->Integral()
           + SBND_int->Integral() + SBND_osc->Integral();

  double nint = SBND_int->Integral();
  double nNC = SBND_NC->Integral();
  double nnumu = SBND_numu->Integral();
  double nosc = SBND_osc->Integral();

  ll.SetFillColor(10);
  ll.AddEntry(SBND_int, Form("Beam #nu_{e} (%3.1f%%)",nint*100./n), "f");
  ll.AddEntry(SBND_NC, Form("NC (%3.1f%%)",nNC*100./n), "f");
  ll.AddEntry(SBND_numu, Form("#nu_{#mu} (%3.1f%%)",nnumu*100./n), "f");
  ll.AddEntry(SBND_osc, Form("Signal #nu_{e} (%3.1f%%)",nosc*100./n), "f");
  // For LSND best-fit values
  ll.AddEntry((TObject*)0, "sin^{2}2#theta_{#mu e} = 0.003","");
  ll.AddEntry((TObject*)0, "#Delta m^{2}_{41} = 1.2 eV^{2}","");

  TLatex* prelim = new TLatex(.9, .96, "SBND Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);

  TCanvas* c1 = new TCanvas("c1");
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();
  c1->SaveAs("output/Stack_nue.root");
  c1->SaveAs("output/Stack_nue.pdf");
  
  TCanvas* c2 = new TCanvas("c2","c2", 700,800);
  c2->Divide(1,2);

  c2->cd(1);
  h->Draw("hist");
  ll.Draw();
  prelim->Draw();

  c2->cd(2);
  h->Draw("nostack hist");
  prelim->Draw();
  c2->SaveAs("output/StackNostack_nue.root");
  c2->SaveAs("output/StackNostack_nue.pdf");
}
