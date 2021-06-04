void build_hist(){
  TFile* beta_file = TFile::Open("beta_file.root");
  auto beta_tree = beta_file->Get<TTree>("beta_tree");

  //create image of all event
  TCanvas* c_all = new TCanvas("c_all", "All");

  TH1F* e_KE_all = new TH1F("e_KE_all", "e_E;Energy [eV];event/bin", 100, 0, 5000000);
  beta_tree->Draw("e_KE>>e_KE_all");

  c_all->SaveAs("all_histogram.png");

  delete c_all;

  //create image of detected event
  TCanvas* c_detected = new TCanvas("c_detected", "Detected");
  
  TH1F* e_KE_detected = new TH1F("e_KE_detected", "e_KE detected;Energy [eV];event/bin", 100, 0, 3000000);
  beta_tree->Draw("e_KE>>e_KE_detected", "e_anihilation_type == 1", "same");
  e_KE_detected->Fit("gaus");

  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.87);
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.20);
  gStyle->SetOptFit();



  c_detected->SaveAs("detected_histogram.png");

  delete c_detected;
  delete e_KE_detected;

  //create image comparing histogram
  TCanvas* c_compare = new TCanvas("c_compare", "Compare");

  e_KE_all->SetAxisRange(0, 3000000, "X");
  beta_tree->Draw("e_KE>>e_KE_all");
  e_KE_all->SetLineColor(kRed);
  beta_tree->Draw("e_KE>>e_KE_detected", "e_anihilation_type == 1", "same");
  e_KE_all->Scale(e_KE_detected->GetBinContent(e_KE_detected->GetMaximumBin())/e_KE_all->GetBinContent(e_KE_all->GetMaximumBin()));

  c_compare->SaveAs("compare_histogram.png");

  delete c_compare;

}