//
//   @file    Belle2Labels.h


#ifndef __BELLE2LABELS_H
#define __BELLE2LABELS_H

#include "Rtypes.h"

#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TPad.h"
#include "TMarker.h"


void BELLE2Label(Double_t x,Double_t y,const char* text,Color_t color = 1) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"Belle II");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
  }
}


void BELLE2LabelOld(Double_t x,Double_t y,bool Preliminary,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"BELLE2");
  if (Preliminary) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+0.115,y,"Preliminary");
  }
}



void BELLE2Version(const char* version,Double_t x,Double_t y,Color_t color) 
{

  if (version) {
    char versionString[100];
    sprintf(versionString,"Version %s",version);
    TLatex l;
    l.SetTextAlign(22); 
    l.SetTextSize(0.04); 
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextColor(color);
    l.DrawLatex(x,y,versionString);
  }
}

#endif // __BELLE2LABELS_H
