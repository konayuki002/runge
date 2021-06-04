void build_hist(){
  TFile* beta_file = TFile::Open("beta_file.root");
  auto beta_tree = beta_file->Get<TTree>("beta_tree");

  //create image of all event
  TCanvas* c_all = new TCanvas("c_all", "All");

  TH1F* e_E_all = new TH1F("e_E_all", "e_E;Energy [eV];event/bin", 100, 0, 5000000);
  beta_tree->Draw("e_E>>e_E_all");

  c_all->SaveAs("all_histogram.png");

  delete c_all;

  //create image of detected event
  TCanvas* c_detected = new TCanvas("c_detected", "Detected");
  
  TH1F* e_E_detected = new TH1F("e_E_detected", "e_E detected;Energy [eV];event/bin", 100, 0, 3000000);
  beta_tree->Draw("e_E>>e_E_detected", "e_anihilation_type == 1", "same");

  c_detected->SaveAs("detected_histogram.png");

  delete c_detected;

  //create image comparing histogram
  TCanvas* c_compare = new TCanvas("c_compare", "Compare");

  e_E_all->SetAxisRange(0, 3000000, "X");
  beta_tree->Draw("e_E>>e_E_all");
  e_E_all->SetLineColor(kRed);
  beta_tree->Draw("e_E>>e_E_detected", "e_anihilation_type == 1", "same");
  e_E_all->Scale(e_E_detected->GetBinContent(e_E_detected->GetMaximumBin())/e_E_all->GetBinContent(e_E_all->GetMaximumBin()));

  c_compare->SaveAs("compare_histogram.png");

  delete c_compare;

}