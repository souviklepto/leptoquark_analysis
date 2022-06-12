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

void crystal_argus_rho() {
  using namespace RooFit;
  SetBelle2Style();
  RooRealVar *deltam = new RooRealVar("deltam","M_{D^{0}}[GeV/c^{2}]", 1.8, 1.9);
  TChain chain("h1");
  chain.Add("/home/souvik/leptoquark_analysis/elec_mode/elec_signalMC/*.root");
  //chain.Add("/home/souvik/leptoquark_analysis/elec_mode/d0bar_signalMC/*.root");
  
  Float_t d0mass,mass,mom,delta,dsmom,dstrf,dz,e1_l,prk1_l,prp1_l,elec_mom,pr_mom,bc0,bc1,pissvd,prsvd,elsvd,pis_k,pis_p ;

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
  chain.SetBranchAddress("pis_p",&pis_p);
  chain.SetBranchAddress("pissvd",&pissvd);
  chain.SetBranchAddress("prsvd",&prsvd);
  chain.SetBranchAddress("elsvd",&elsvd);
  chain.SetBranchAddress("bc0",&bc0);
  chain.SetBranchAddress("bc1",&bc1);
  
  RooDataSet* data=new RooDataSet("data","",RooArgSet(*deltam));
  for(int i=0;i< chain.GetEntries();i++) {
    chain.GetEntry(i);
    double q=mass-d0mass-0.139 ;                                                    
    if(d0mass>1.8 && d0mass<1.9 && delta>0.144 && delta<0.147 && e1_l>0.9 && prk1_l>0.6 && prp1_l>0.6 && q<0.02 && mom>2.5 && pis_p>0.6 && pissvd>2 && prsvd>2 && elsvd>2 && elec_mom>0.7 && pr_mom>1.7 && bc0==100)
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
  RooRealVar *sls1=new RooRealVar("sls1","Ratio of Sigma",0.,15);
  RooRealVar *srs1=new RooRealVar("srs1","Ratio of Sigma",0.,30);
  RooRealVar *a2a=new RooRealVar("a2a", " Area2/Area1",.01,1);
  RooRealVar *delm1=new RooRealVar("delm1", "Difference of Mean",0);
  RooFormulaVar *mean2= new RooFormulaVar("mean3","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1,*delm1));
  RooFormulaVar *sigmaL= new RooFormulaVar("sigmaL", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1,*sls1));
  RooFormulaVar *sigmaR= new RooFormulaVar("sigmaR", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1,*srs1));
  RooBifurGauss *gauss3_1= new RooBifurGauss("gauss3_1", "3nd gaussian PDF", *deltam, *mean2, *sigmaL,*sigmaR);
  RooAddPdf *dG_1= new RooAddPdf("dG_1", " 1st GAuss + 3rd Gauss", RooArgList(*gauss_cb,*gauss), RooArgList(*a2a));
  
  RooCBShape *gauss1_1 = new RooCBShape("gauss1_1","Crystal Ball shape", *deltam, *mean, *sigma, *alpha, *n);    

  RooRealVar *argpar = new RooRealVar("argpar","argus shape parameter",-80,1);
  RooRealVar *chi = new RooRealVar("chi","argus shape parameter",5.29);
  RooArgusBG *argus = new RooArgusBG("argus","Argus PDF",*deltam,*chi,*argpar);

  RooRealVar * slope1= new RooRealVar("slope1", "Slope of Polynomial", -10, 10);
  RooRealVar * slope2= new RooRealVar("slope2", " Slope of Cheb2", -1, 1);

  RooChebychev *chebpol = new RooChebychev("chebpol","Chebshev Polynomial ", *deltam, RooArgList(*slope1,*slope2));

  //RooRealVar *a2a=new RooRealVar("a2a", " Area2/Area1",0.01,1);

  RooRealVar *sig = new RooRealVar("sig", "",8000,200000);
  RooRealVar *BKG = new RooRealVar("BKG", "",3500,300000);

  //RooAddPdf *dG_1= new RooAddPdf("dG_1", " 2nd Gauss + Argus", RooArgList(*gauss2_1,*gauss3_1), RooArgList(*a2a));

  RooAddPdf *depdf = new RooAddPdf ("depdf", "Two Gaussian +  ",RooArgList(*dG_1,*chebpol), RooArgList(*sig,*BKG));

  TCanvas *c1 =new TCanvas ("c1", "D^{0} Mass", 700, 600);
  c1->Divide(2,2);
  RooFitResult *fitresult = depdf->fitTo(*data, Extended(true),Minos(true));
  RooPlot* deltaEplot = deltam->frame(Bins(50));
  data->plotOn(deltaEplot);
  depdf->paramOn(deltaEplot);
  depdf->plotOn(deltaEplot);
  //depdf->plotOn(deltaEplot,Components(RooArgSet(*chebpol)),LineColor(kGreen),LineStyle(kDashed));                           
  //depdf->plotOn(deltaEplot,Components(RooArgSet(*argus)),LineColor(kCyan),LineStyle(kDashed));
  depdf->plotOn(deltaEplot,Components(RooArgSet(*chebpol)),LineColor(kRed),LineStyle(kDashed));                               
  depdf->plotOn(deltaEplot,Components(RooArgSet(*gauss3_1)),LineColor(kMagenta),LineStyle(kDashed));                                       
  depdf->plotOn(deltaEplot,Components(RooArgSet(*gauss)),LineColor(kGreen),LineStyle(kDashed));
  depdf->plotOn(deltaEplot,Components(RooArgSet(*dG_1)),LineColor(kBlue),LineStyle(kDashed));                                     
  deltaEplot->SetTitle("");
  deltaEplot->GetYaxis()->CenterTitle();
  deltaEplot->GetXaxis()->CenterTitle();
  deltaEplot->Draw();
  
}
