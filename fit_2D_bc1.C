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

void fit_2D_bc1()
{
  using namespace RooFit;
  SetBelle2Style();
  /*******************Fit Variables***********************************/

  RooRealVar *deltam = new RooRealVar("deltam","", 0.14, 0.152);
  RooRealVar *md0 = new RooRealVar("md0","", 1.8, 1.9);
  
  /*******************Input root file**********************************/
  
  TChain chain("h1");
  chain.Add("/home/souvik/leptoquark_analysis/elec_mode/d0barempp_mode/evtgen_exp_*.root");
  Float_t d0mass,mass,mom,delta,dsmom,dstrf,dz,e1_l,prk1_l,prp1_l,elec_mom,pr_mom,pissvd,elsvd,prsvd,pis_k,bc3,bc1,pis_p ;

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
  chain.SetBranchAddress("bc3",&bc3);
  chain.SetBranchAddress("bc1",&bc1);
    
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(*deltam,*md0));
  for(int i=0;i< chain.GetEntries();i++) {
    chain.GetEntry(i);
    double q=mass-d0mass-0.139 ;
    deltam->setVal(delta);
    md0->setVal(d0mass);
    if(d0mass>1.8 && d0mass<1.9 && delta>0.14 && delta<0.152 && e1_l>0.9 && prk1_l>0.6 && prp1_l>0.6 && q<0.02 && mom>2.5 && pis_p>0.6 && pissvd>2 && prsvd>2 && elsvd>2 && elec_mom>0.7 && pr_mom>1.7 && bc1==100)
      data->add(RooArgSet(*deltam,*md0));
  }
  
  /*****************************Fit***********************/
  RooRealVar *mean1_dm = new RooRealVar("mean1_dm", "MEAN of 1st gaussian",0.145,0.14,0.152);
  RooRealVar *sigma1_dm = new RooRealVar("sigma1_dm", "Sigma of 1st gaussian",0.0002618,0,0.0089);
  RooGaussian *gauss_dm= new RooGaussian("gauss_dm", "1st gaussian PDF", *deltam, *mean1_dm, *sigma1_dm);        
  RooRealVar *delm1_dm=new RooRealVar("delm1_dm", "Difference of Mean",-0.00018805,-0.01,0.01);
  RooFormulaVar *mean2_dm= new RooFormulaVar("mean2_dm","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1_dm,*delm1_dm));
  RooRealVar *mean3_dm = new RooRealVar("mean3_dm", "MEAN of 1st gaussian",0.145,0.14,0.152);
  RooRealVar *sigma2_dm = new RooRealVar("sigma2_dm", "Sigma of 2nd gaussian",0.00298,0,0.0089);
  RooRealVar *sigma3_dm = new RooRealVar("sigma3_dm", "Sigma of 2nd gaussian",0.000504,0,0.0089);
  RooRealVar *sigma4_dm = new RooRealVar("sigma4_dm", "Sigma of 2nd gaussian",0.001167,0,0.0089);
  RooRealVar *mean4_dm = new RooRealVar("mean4_dm", "MEAN of 1st gaussian",0.145,0.14,0.152);
  RooGaussian *gauss2_dm= new RooGaussian("gauss2_dm", "1st gaussian PDF", *deltam, *mean2_dm, *sigma2_dm);
  RooGaussian *gauss3_dm= new RooGaussian("gauss3_dm", "3rd gaussian PDF", *deltam, *mean3_dm, *sigma3_dm);
  RooGaussian *gauss4_dm= new RooGaussian("gauss4_dm", "4th gaussian PDF", *deltam, *mean4_dm, *sigma4_dm);
  
  RooRealVar *a2a_dm1=new RooRealVar("a2a_dm1", " Area2/Area1",0.9351,0.1,1);
  RooRealVar *a2a_dm=new RooRealVar("a2a_dm", " Area2/Area1",0.606,0.1,1);
  RooRealVar *a2a_dm2=new RooRealVar("a2a_dm2", " Area2/Area1",0.9244,0.01,1);
  RooAddPdf *dG_dm3= new RooAddPdf("dG_dm3", " 1st GAuss + 2nd Gauss", RooArgList(*gauss_dm,*gauss2_dm), RooArgList(*a2a_dm1));
  RooAddPdf *dG_dm2= new RooAddPdf("dG_dm2", " 1st GAuss + 2nd Gauss", RooArgList(*dG_dm3,*gauss3_dm), RooArgList(*a2a_dm));  
  RooAddPdf *dG_dm4= new RooAddPdf("dG_dm4", " 1st GAuss + 2nd Gauss", RooArgList(*dG_dm2,*gauss4_dm), RooArgList(*a2a_dm2));

  //------------------BKG PDFs--------------------                                                                                                                
  RooRealVar th("th","a",0.140);
  RooRealVar a("a","a",-436,-5000,5000);//-436,-5000,5000);
  RooRealVar b("b","",-4102,-2000000,2000000);//-41023,-200000,200000);
  RooRealVar c("c","",0);//177,-100000,100000);
  RooRealVar d("d","",0);//-5,5)                                                                                                  
  RooGenericPdf *t_bkg = new RooGenericPdf("t_bkg",
                                           "(deltam>=th)*pow((deltam-th),4)*exp(a*(deltam-th)+ b*pow((deltam-th),2))",
                                           RooArgSet(*deltam,th,a,b));
  RooRealVar * slope4= new RooRealVar("slope4", " Slope of Cheb2", -1,1);
  RooRealVar * slope5= new RooRealVar("slope5", " Slope of Cheb2", -1,1);
  RooChebychev *chebpol1 = new RooChebychev("chebpol1","Chebshev Polynomial ", *deltam, RooArgList(*slope4,*slope5));
  RooRealVar *a2a_dm3=new RooRealVar("a2a_dm3", " Area2/Area1",0.01,1);
  RooAddPdf *dG_dmb= new RooAddPdf("dG_dmb", "bkg1 + bkg2", RooArgList(*t_bkg,*chebpol1), RooArgList(*a2a_dm3));
    
  // M_D0 fit
  RooRealVar *mean1_md = new RooRealVar("mean1_md", "MEAN of 1st gaussian",1.864693,1.8,1.9);
  RooRealVar *sigma1_md = new RooRealVar("sigma1_md", "Sigma of 1st gaussian",0.003236,0,0.0089);
  RooGaussian *gauss_md= new RooGaussian("gauss_md", "1st gaussian PDF", *md0, *mean1_md, *sigma1_md);
  RooRealVar *mean_md = new RooRealVar("mean_md", "MEAN of 1st gaussian",1.8,1.9);
  RooRealVar *sigma_md = new RooRealVar("sigma_md", "Sigma of 1st gaussian",0.00716,0,0.05);
  RooGaussian *gauss_mds= new RooGaussian("gauss_mds", "1st gaussian PDF", *md0, *mean1_md, *sigma_md);
  RooRealVar *sls1_md=new RooRealVar("sls1_md","Ratio of Sigma",0.1,20);
  RooRealVar *srs1_md=new RooRealVar("srs1_md","Ratio of Sigma",0.95,0.1,10);
  RooRealVar *a2a_md=new RooRealVar("a2a_md", " Area2/Area1",0.1,1);
  RooRealVar *delm1_md=new RooRealVar("delm1_md", "Difference of Mean",0.00484,-0.01,0.01);
  RooFormulaVar *mean2_md= new RooFormulaVar("mean2_md","Mean of 2nd gaussian", "@0+@1",RooArgList(*mean1_md,*delm1_md));
  RooFormulaVar *sigmaL_md= new RooFormulaVar("sigmaL_md", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1_md,*sls1_md));
  RooFormulaVar *sigmaR_md= new RooFormulaVar("sigmaR_md", "Sigma of 2nd gaussian", " @0*@1", RooArgList(*sigma1_md,*srs1_md));
  RooBifurGauss *gauss3_md= new RooBifurGauss("gauss3_md", "3nd gaussian PDF", *md0, *mean2_md, *sigmaL_md,*sigmaR_md);
  RooAddPdf *dG_md1= new RooAddPdf("dG_md1", " 1st GAuss + 3rd Gauss", RooArgList(*gauss3_md,*gauss_md), RooArgList(*a2a_md));
  RooRealVar *a2a_md1=new RooRealVar("a2a_md1", " Area2/Area1",0.761,0.1,1);
  RooAddPdf *dG_md2= new RooAddPdf("dG_md2", " 2nd Gauss + Argus", RooArgList(*dG_md1,*gauss_mds), RooArgList(*a2a_md1));
  
  RooRealVar * slope1= new RooRealVar("slope1", "Slope of Polynomial",-0.3197,-1,1);
  RooRealVar * slope2= new RooRealVar("slope2", " Slope of Cheb2",-0.4897,-1,1);
  RooRealVar * slope3= new RooRealVar("slope3", "Slope of Polynomial", -1,1);
  RooChebychev *chebpol = new RooChebychev("chebpol","Chebshev Polynomial ", *md0, RooArgList(*slope1,*slope2));
  
  //RooProdPdf *pdf_signal = new RooProdPdf("pdf_signal","pdf_signal",RooArgList(*t_bkg,*chebpol1));  //2D signal PDF                                                                                    
  RooProdPdf *pdf_signal_signal = new RooProdPdf("pdf_signal_signal","pdf_signal_signal",RooArgList(*dG_dm4,*dG_md2));
  RooProdPdf *pdf_signal = new RooProdPdf("pdf_signal","pdf_signal",RooArgList(*t_bkg,*chebpol));
  
  RooRealVar *sig = new RooRealVar("sig", "",8000,200000);
  RooRealVar *BKG = new RooRealVar("BKG", "",200,8000);
  RooAddPdf *total = new RooAddPdf("total","total",RooArgList(*pdf_signal_signal,*pdf_signal),RooArgList(*sig,*BKG));
    
  total->fitTo(*data);

  deltam->setRange("signal1",0.144,0.147);
  md0->setRange("signal2",1.85,1.875);
  
  RooPlot *xframe = data->plotOn(deltam->frame(100),CutRange("signal2"));
  //RooPlot *xframe = data->plotOn(deltam.frame(50));                                                                                                                                                      
  total->plotOn(xframe,ProjectionRange("signal2"));

  Double_t chisq = xframe->chiSquare();

  total->paramOn(xframe,data);

  RooHist* hpull = xframe->pullHist() ;

  total->plotOn(xframe, Components(RooArgList(*pdf_signal_signal)),LineStyle(kDashed),LineColor(kGreen),ProjectionRange("signal2"));
  total->plotOn(xframe, Components(RooArgList(*pdf_signal)),LineStyle(kDashed),LineColor(kRed),ProjectionRange("signal2"));

  RooPlot* frame3 = deltam->frame(Title("")) ;
  frame3->addPlotable(hpull,"P") ;
  TCanvas* c1 = new TCanvas() ;
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->Draw();             // Draw the upper pad: pad1                                                                                                                                                   
  pad1->cd();
  xframe->GetYaxis()->CenterTitle();
  xframe->SetTitle("");
  xframe->Draw() ;

  TLatex s;
     s.SetNDC();
     s.SetTextFont(42);
     s.SetTextColor(2);
     s.DrawLatex(0.3,0.5,"#chi^{2}/ndf = 1.33");
  
  c1->cd();          // Go back to the main canvas before defining pad2                                                                                                                                    
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->Draw();
  pad2->cd();
  frame3->GetXaxis()->SetLabelSize(0.15);
  frame3->GetXaxis()->SetTitleSize(0.15);
  frame3->GetXaxis()->SetTitle("M_{D^{*}}-M_{D^{0}}[GeV/c^{2}]");
  frame3->GetXaxis()->CenterTitle();
  frame3->GetYaxis()->SetNdivisions(105);
  frame3->GetYaxis()->SetTickLength(0.075);
  frame3->GetYaxis()->SetTitleSize(0.15);
  frame3->GetYaxis()->SetLabelSize(0.15);
  frame3->GetYaxis()->SetTitleOffset(0.3);
  frame3->GetYaxis()->SetTitle("Pull");
  frame3->GetYaxis()->CenterTitle();
  frame3->SetTitle("");
  frame3->Draw() ;
  cout << "chi2 deltam=" << chisq << endl;

  //Plotting MD projection                                                                                                                                                                                
  RooPlot *yframe = data->plotOn(md0->frame(100),CutRange("signal1"));
  total->plotOn(yframe,ProjectionRange("signal1"));

  Double_t chisq_md0 = yframe->chiSquare();
  RooHist* hpull1 = yframe->pullHist() ;

  total->plotOn(yframe, Components(RooArgList(*pdf_signal_signal)),LineStyle(kDashed),LineColor(kGreen),ProjectionRange("signal1"));
  total->plotOn(yframe, Components(RooArgList(*pdf_signal)),LineStyle(kDashed),LineColor(kRed),ProjectionRange("signal2"));
    
  RooPlot* frame5 = md0->frame(Title("")) ;
  frame5->addPlotable(hpull1,"P");
  TCanvas* c2 = new TCanvas() ;
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->Draw();                                                                                                                         
  pad11->cd();
  yframe->GetYaxis()->CenterTitle();
  yframe->SetTitle("");
  yframe->Draw() ;

   TLatex l;
     l.SetNDC();
     l.SetTextFont(42);
     l.SetTextColor(2);
     l.DrawLatex(0.3,0.5,"#chi^{2}/ndf = 1.36");
  
  c2->cd();                                                                                                                                   
  TPad *pad22 = new TPad("pad22", "pad22", 0, 0.05, 1, 0.3);
  pad22->Draw();
  pad22->cd();
  frame5->GetXaxis()->SetLabelSize(0.15);
  frame5->GetXaxis()->SetTitleSize(0.15);
  frame5->GetXaxis()->CenterTitle();
  frame5->GetXaxis()->SetTitle("M_{D^{0}}[GeV/c^{2}]");
  frame5->GetYaxis()->SetNdivisions(105);
  frame5->GetYaxis()->SetTickLength(0.075);
  frame5->GetYaxis()->SetTitleSize(0.15);
  frame5->GetYaxis()->SetLabelSize(0.15);
  frame5->GetYaxis()->SetTitleOffset(0.3);
  frame5->GetYaxis()->SetTitle("Pull");
  frame5->GetYaxis()->CenterTitle();
  frame5->SetTitle("");
  frame5->Draw() ;
  cout << "chi2 md0=" << chisq_md0 << endl;
  
}
