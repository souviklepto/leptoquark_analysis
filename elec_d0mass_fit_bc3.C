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

void elec_d0mass_fit_bc3() {
  using namespace RooFit;
  SetBelle2Style();
  RooRealVar *deltam = new RooRealVar("deltam","M_{D^{0}}[GeV/c^{2}]", 1.8, 1.9);
  TChain chain("h1");
  chain.Add("/home/souvik/leptoquark_analysis/elec_mode/d0barpmep_mode/evtgen_exp_*.root");
  Float_t d0mass,mass,mom,delta,dsmom,dstrf,dz,e1_l,prk1_l,prp1_l,elec_mom,pr_mom,bc2,bc3,pissvd,prsvd,elsvd,pis_k,pis_p ;
  chain.SetBranchAddress("d0mass",&d0mass);
  chain.SetBranchAddress("mass",&mass);
  chain.SetBranchAddress("mom",&mom);
  chain.SetBranchAddress("deltam",&delta);
  chain.SetBranchAddress("dz",&dz);
  chain.SetBranchAddress("dstrf",&dstrf);
  chain.SetBranchAddress("e1_l",&e1_l);
  chain.SetBranchAddress("prk1_l",&prk1_l);
  chain.SetBranchAddress("prp1_l",&prp1_l);
  chain.SetBranchAddress("elec_mom",&elec_mom);
  chain.SetBranchAddress("pr_mom",&pr_mom);
  chain.SetBranchAddress("pissvd",&pissvd);
  chain.SetBranchAddress("elsvd",&elsvd);
  chain.SetBranchAddress("prsvd",&prsvd);
  chain.SetBranchAddress("pis_p",&pis_p);
  chain.SetBranchAddress("bc2",&bc2);
  chain.SetBranchAddress("bc3",&bc3);
  
  RooDataSet* data=new RooDataSet("data","",RooArgSet(*deltam));
  for(int i=0;i< chain.GetEntries();i++) {
    chain.GetEntry(i);
    double q=mass-d0mass-0.139 ;                                                    
    if(d0mass>1.8 && d0mass<1.9 && delta>0.144 && delta<0.147 && pissvd>2 && elsvd>2 && prsvd>2 && e1_l>0.9 && prk1_l>0.6 && prp1_l>0.6 && pis_p>0.6 && mom>2.5 && pr_mom>1.7 && elec_mom>0.7 && bc3==100)
     {deltam->setVal(d0mass);
      data->add(RooArgSet(*deltam));
    }
  }

  RooRealVar *mean = new RooRealVar("mean", "MEAN of 1st gaussian",1.8,1.9);
  RooRealVar *sigma = new RooRealVar("sigma", "Sigma of 1st gaussian",0,0.03);
  RooRealVar *alpha = new RooRealVar("alpha", "alpha",3,0,5);
  RooRealVar *n = new RooRealVar("n", "Polynom degree",4,1,15);
  RooCBShape *gauss_cb= new RooCBShape("gauss_cb", "1st gaussian PDF", *deltam, *mean, *sigma, *alpha, *n);

  RooRealVar *mean1 = new RooRealVar("mean1", "MEAN of 1st gaussian",1.8,1.9);
  RooRealVar *sigma1 = new RooRealVar("sigma1", "Sigma of 1st gaussian",0,0.005);
  RooGaussian *gauss= new RooGaussian("gauss", "1st gaussian PDF", *deltam, *mean1, *sigma1);
  RooRealVar *mean_s = new RooRealVar("mean_s", "MEAN of 1st gaussian",1.8,1.9);
  RooRealVar *sigma_s = new RooRealVar("sigma_s", "Sigma of 1st gaussian",0,0.01);
  RooGaussian *gauss_s= new RooGaussian("gauss_s", "1st gaussian PDF", *deltam, *mean1, *sigma_s);
  RooRealVar *sls1=new RooRealVar("sls1","Ratio of Sigma",0.,15);
  RooRealVar *srs1=new RooRealVar("srs1","Ratio of Sigma",0.,30);
  RooRealVar *a2a=new RooRealVar("a2a", " Area2/Area1",.01,1);
  RooRealVar *delm1=new RooRealVar("delm1", "Difference of Mean",-0.01,0.01);
  RooFormulaVar *mean2= new RooFormulaVar("mean3","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1,*delm1));
  RooFormulaVar *sigmaL= new RooFormulaVar("sigmaL", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1,*sls1));
  RooFormulaVar *sigmaR= new RooFormulaVar("sigmaR", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1,*srs1));
  RooBifurGauss *gauss3_1= new RooBifurGauss("gauss3_1", "3nd gaussian PDF", *deltam, *mean2, *sigmaL,*sigmaR);
  RooAddPdf *dG_1= new RooAddPdf("dG_1", " 1st GAuss + 3rd Gauss", RooArgList(*gauss3_1,*gauss), RooArgList(*a2a));
  
  RooCBShape *gauss1_1 = new RooCBShape("gauss1_1","Crystal Ball shape", *deltam, *mean, *sigma, *alpha, *n);    

  RooRealVar * slope1= new RooRealVar("slope1", "Slope of Polynomial", -10, 10);
  RooRealVar * slope2= new RooRealVar("slope2", " Slope of Cheb2", -1, 1);
  RooChebychev *chebpol = new RooChebychev("chebpol","Chebshev Polynomial ", *deltam, RooArgList(*slope1,*slope2));

  RooRealVar *a2a_1=new RooRealVar("a2a_1", " Area2/Area1",0.01,1);

  RooRealVar *sig = new RooRealVar("sig", "",800,200000);
  RooRealVar *BKG = new RooRealVar("BKG", "",350,300000);

  RooAddPdf *dG_2= new RooAddPdf("dG_2", " 2nd Gauss + Argus", RooArgList(*dG_1,*gauss_s), RooArgList(*a2a_1));
  RooAddPdf *depdf = new RooAddPdf ("depdf", "Two Gaussian +  ",RooArgList(*dG_2,*chebpol), RooArgList(*sig,*BKG));

  TCanvas *c1 =new TCanvas ("c1", "D^{0} Mass", 700, 600);
  c1->Divide(2,2);
  RooFitResult *fitresult = depdf->fitTo(*data, Extended(true),Minos(true));
  RooPlot* frame1 = deltam->frame(Bins(100));
  data->plotOn(frame1);
  depdf->paramOn(frame1);
  depdf->plotOn(frame1);
  double chi2_ndf1 =  frame1->chiSquare();
  cout << "chi2_ndf1 in data1 = " << chi2_ndf1 << endl;
  //depdf->plotOn(frame1);                                                                                                                                                                                 
  RooHist* hpull = frame1->pullHist();
  hpull->SetFillColor(kAzure-8);
  RooPlot* frame2 = deltam->frame(Title(" ")) ;
  frame2->addPlotable(hpull,"P") ;
  RooArgSet* params = depdf->getParameters(*deltam) ;
  depdf->plotOn(frame1,Components(RooArgSet(*chebpol)),LineColor(kRed),LineStyle(kDashed));                               
  depdf->plotOn(frame1,Components(RooArgSet(*gauss3_1)),LineColor(kViolet+6),LineStyle(kDashed));                                       
  depdf->plotOn(frame1,Components(RooArgSet(*gauss)),LineColor(kGreen),LineStyle(kDashed));
  depdf->plotOn(frame1,Components(RooArgSet(*gauss_s)),LineColor(kCyan+2),LineStyle(kDashed));
  //depdf->plotOn(frame1,Components(RooArgSet(*dG_1)),LineColor(kBlue),LineStyle(kDashed));                                     
  frame1->SetTitle("");
  frame1->GetYaxis()->CenterTitle();
  frame1->GetXaxis()->CenterTitle();
  TPad *Pad1 = new TPad("Pad1", " ",0.0,0.25,1.0,1.0);
  Pad1->Draw();
  TPad *Pad2 = new TPad("Pad2", " ",0.0,0.0,1.0,0.25);
  Pad2->Draw();
  Pad1->cd();
  Pad1->SetTopMargin(0.08);
  Pad1->SetRightMargin(0.05);
  Pad1->SetLeftMargin(0.1);
  Pad1->SetBottomMargin(0.025);
  frame1->GetXaxis()->SetLabelSize(0.03);
  frame1->GetXaxis()->SetTitleSize(0.0);
  frame1->GetXaxis()->SetLabelSize(0.0);
  frame1->GetYaxis()->SetLabelSize(0.03);
  frame1->GetYaxis()->SetTitleSize(0.05);
  frame1->GetYaxis()->SetTitleOffset(0.9);
  //frame1->GetYaxis()->SetMaxDigits(3);                                                                                                                                                                 
  frame1->SetStats(0);
  frame1->Draw();

  TLatex s;
  s.SetNDC();
  s.SetTextFont(42);
  s.SetTextColor(1);
  s.DrawLatex(0.3,0.5,"#chi^{2}/ndf = 1.05");

  Pad2->cd();
  Pad2->SetTopMargin(0.03);
  Pad2->SetRightMargin(0.05);
  Pad2->SetLeftMargin(0.1);
  Pad2->SetBottomMargin(0.4);

  frame2->GetXaxis()->SetNdivisions(505);
  frame2->GetXaxis()->SetTickLength(0.1);
  frame2->GetXaxis()->SetTitleSize(0.15);
  frame2->GetXaxis()->SetLabelSize(0.15);
  frame2->GetXaxis()->SetTitleOffset(0.95);

   frame2->GetXaxis()->SetTitle("M_{D^{0}}[GeV/c^{2}]");
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetXaxis()->CenterTitle();
  frame2->GetYaxis()->SetNdivisions(105);
  frame2->GetYaxis()->SetTickLength(0.075);
  frame2->GetYaxis()->SetTitleSize(0.15);
  frame2->GetYaxis()->SetLabelSize(0.15);
  frame2->GetYaxis()->SetTitleOffset(0.3);
  //frame2->GetYaxis()->SetMaxDigits(3);                                                                                                                                                                 
  frame2->Draw();

  
  
}
