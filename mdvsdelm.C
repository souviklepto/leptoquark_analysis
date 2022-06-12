{

  gStyle->SetOptStat("emr");

  TChain *ch = new TChain("h1");
  ch->Add("/home/souvik/new_belle_analysis/elec_mode/elec_signalMC/evtgen_*.root");
  TTree *tr = ch;
  
  TH2D *hist   = new TH2D("hist"  ,"pr_mom vs mom ; pr_mom ; mom",75, 0, 5, 75, 0, 5);

  //TCut m_p = "mom > 2.5";
  
  hist->SetMarkerStyle(20);
  hist->GetYaxis()->SetTitleOffset(1.6);
  //hist->Scale(0.10);
  TCanvas *c1 = new TCanvas("c1","",500,500);
  //tr->Draw("e1_l : prp1_l >> hist",m_p);
  tr->Draw("pr_mom : mom >> hist");
  
  //hist->GetCovariance(Double_t axis=1 , Double_t axis=2);

  auto cov = hist->GetCovariance(2,1);
  std::cout << "covariance: " << cov << std::endl;
  
}
