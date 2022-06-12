#include "Belle2Style.h"
#include "Belle2Utils.h"
#include "Belle2Labels.h"

void SetBelle2Style() {
  static TStyle* belle2Style = 0;
  cout << "\nApplying BELLE2 style settings...\n" << endl ;
  if (!belle2Style)
    belle2Style = Belle2Style();
  gROOT->SetStyle("BELLE2");
  gROOT->ForceStyle();
}

void deltam_tm()
{
  SetBelle2Style();
  TChain *t1 = new TChain();
  t1->Add("/home/souvik/leptoquark_analysis/elec_mode/d0bar_signalMC/evtgen_*.root/h1");
  Float_t mult_ds, dstrf, d0f,d0mass,mass,mom,deltam,elf,prf,elec_mom,pr_mom,elec_mas,pr_mass,mult_d0,bc0,e1_l,prk1_l,prp1_l,pis_k,pissvd,prsvd,elsvd,bc1 ;
  
  t1->SetBranchAddress("mult_ds",&mult_ds);                                                        
  t1->SetBranchAddress("e1_l",&e1_l);
  t1->SetBranchAddress("prk1_l",&prk1_l);
  t1->SetBranchAddress("prp1_l",&prp1_l);
  t1->SetBranchAddress("pis_k",&pis_k);                                                              
  t1->SetBranchAddress("deltam",&deltam);                                                          
  t1->SetBranchAddress("elec_mas",&elec_mas);                                                      
  t1->SetBranchAddress("pr_mass",&pr_mass);                                                      
  t1->SetBranchAddress("elec_mom",&elec_mom);                                                        
  t1->SetBranchAddress("pr_mom",&pr_mom);                                                        
  t1->SetBranchAddress("mom",&mom);                                                                
  t1->SetBranchAddress("dstrf",&dstrf);                                                              
  t1->SetBranchAddress("d0mass",&d0mass);                                                          
  t1->SetBranchAddress("mass",&mass);
  t1->SetBranchAddress("pissvd",&pissvd);
  t1->SetBranchAddress("prsvd",&prsvd);
  t1->SetBranchAddress("elsvd",&elsvd);
  t1->SetBranchAddress("bc1",&bc1);
  
  TH1F *hist3 = new TH1F("hist3","#Delta M",80, 0.14,0.16);
  for (Int_t i=0; i < t1->GetEntries(); ++i)
    {
      t1->GetEntry(i);
      double q = mass-d0mass-0.139 ;
      {
	if(d0mass>1.8 && d0mass<1.9 && mom>2.5 && q<0.02 && e1_l>0.9 && prk1_l>0.9 && prp1_l>0.9 && pis_k<0.4 && pissvd>2 && prsvd>2 && elsvd>2 && bc1==100)
	hist3->Fill(deltam);
      }
    }

  hist3->GetYaxis()->SetRangeUser(0,50000);
  hist3->GetXaxis()->SetTitle("#DeltaM[GeV/c^{2}]");
  hist3->GetYaxis()->SetTitle("Events");
  hist3->GetYaxis()->CenterTitle();
  hist3->GetXaxis()->CenterTitle();
  hist3->SetLineColor(kBlue);
  hist3->SetLineWidth(2);
  hist3->SetFillColor(kBlue);
  hist3->SetFillStyle(3040);
  hist3->Draw();
  
  TH1F *hist4 = new TH1F("hist4","#Delta M",80, 0.14,0.16);
  for (Int_t i=0; i < t1->GetEntries(); ++i)
    {
      t1->GetEntry(i);
      double q = mass-d0mass-0.139 ;
      {
	if(dstrf==1 && d0mass>1.8 && d0mass<1.9 && mom>2.5 && q<0.02 && e1_l>0.9 && prk1_l>0.9 && prp1_l>0.9 && pis_k<0.4 && pissvd>2 && prsvd>2 && elsvd>2 && bc1==100)
        hist4->Fill(deltam);
      }
    }

  hist4->GetYaxis()->SetRangeUser(0,50000);
  hist4->GetXaxis()->SetTitle("#DeltaM[GeV/c^{2}]");
  hist4->GetYaxis()->SetTitle("Events");
  hist4->GetYaxis()->CenterTitle();
  hist4->GetXaxis()->CenterTitle();
  hist4->SetLineColor(kGreen);
  hist4->SetLineWidth(2);
  hist4->SetFillColor(kGreen);
  hist4->SetFillStyle(3004);
  hist4->Draw("SAMES");
  

  TH1F *hist5 = new TH1F("hist5","#Delta M",80, 0.14,0.16);
  for (Int_t i=0; i < t1->GetEntries(); ++i)
    {
      t1->GetEntry(i);
      double q = mass-d0mass-0.139 ;
      {
        if(dstrf!=1 && d0mass>1.8 && d0mass<1.9 && mom>2.5 && q<0.02 && e1_l>0.9 && prk1_l>0.9 && prp1_l>0.9 && pis_k<0.4 && pissvd>2 && prsvd>2 && elsvd>2 && bc1==100)
	  hist5->Fill(deltam);
      }
    }
  
  hist5->GetYaxis()->SetRangeUser(0,50000);
  hist5->GetXaxis()->SetTitle("#DeltaM[GeV/c^{2}]");
  hist5->GetYaxis()->SetTitle("Events");
  hist5->GetYaxis()->CenterTitle();
  hist5->GetXaxis()->CenterTitle();
  hist5->SetLineColor(kRed);
  hist5->SetLineWidth(2);
  hist5->SetFillColor(kRed);
  hist5->SetFillStyle(3003);
  hist5->Draw("SAMES");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(hist3,"tot","l");
  legend->AddEntry(hist4,"sig","l");
  legend->AddEntry(hist5,"scf","l");
  legend->Draw();
  
}

