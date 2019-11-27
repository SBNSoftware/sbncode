{
  // event_numbers.C has to be executed first!
  // Replace Icarus with the appropriate experiment

  TFile* f = new TFile("output/output_nue_Icarus.root");
  
  TH1D* Icarus_NC = (TH1D*) f->Get("Icarus_NC");
  TH1D* Icarus_numu = (TH1D*) f->Get("Icarus_numu");
  TH1D* Icarus_int = (TH1D*) f->Get("Icarus_int");
  TH1D* Icarus_osc = (TH1D*) f->Get("Icarus_osc");

  Icarus_NC->SetFillColorAlpha(Icarus_NC->GetLineColor(), 0.6);
  Icarus_numu->SetFillColorAlpha(Icarus_numu->GetLineColor(), 0.6);
  Icarus_int->SetFillColorAlpha(Icarus_int->GetLineColor(), 0.6);
  Icarus_osc->SetFillColorAlpha(Icarus_osc->GetLineColor(), 0.6);

  //Change POT if needed
  THStack* h = new THStack("h",";Reconstructed energy (GeV);Events / (6.6#times10^{20} POT)");

  h->Add(Icarus_int);
  h->Add(Icarus_NC);
  h->Add(Icarus_numu);
  h->Add(Icarus_osc);

  TLegend ll(0.65,0.5,0.84,0.9);

  double n = Icarus_NC->Integral() + Icarus_numu->Integral()
           + Icarus_int->Integral() + Icarus_osc->Integral();

  double nint = Icarus_int->Integral();
  double nNC = Icarus_NC->Integral();
  double nnumu = Icarus_numu->Integral();
  double nosc = Icarus_osc->Integral();

  ll.SetFillColor(10);
  ll.AddEntry(Icarus_int, Form("Beam #nu_{e} (%3.1f%%)",nint*100./n), "f");
  ll.AddEntry(Icarus_NC, Form("NC (%3.1f%%)",nNC*100./n), "f");
  ll.AddEntry(Icarus_numu, Form("#nu_{#mu} (%3.1f%%)",nnumu*100./n), "f");
  ll.AddEntry(Icarus_osc, Form("Signal #nu_{e} (%3.1f%%)",nosc*100./n), "f");

  TLatex* prelim = new TLatex(.9, .96, "Icarus Simulation");
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
