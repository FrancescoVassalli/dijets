#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TChain.h"

/**
macro for making the histograms based on fat ratio

**/


void FatHist(int filecount){
	TCanvas *tc = new TCanvas();
 	bool lowpT = true;
 	float stateMax;
  	std::string outfile;
  	if(lowpT){
   		stateMax=3.5;
    	outfile = " 101.pdf";
  	}
  	else{
    	stateMax = 4.5;
    	outfile = " 200.pdf";
  	}
  	TChain *dijet_tree = new TChain("tree100");

  	std::string infile;
  	if(lowpT) 
  		infile= "Rjet11";
  	else
  		infile= "Rjet20";
  	std::string fileappend = ".root";
  	int printcount=0;
  	std::string filename = infile+std::to_string(printcount)+fileappend;
  	while(printcount<=filecount){
  		dijet_tree->Add(filename.c_str());
  		printcount++;
  		filename=infile+std::to_string(printcount)+fileappend;
  	}
  	gStyle->SetOptStat(0);
  	dijet_tree->Add(infile.c_str());

  	float bins[] = {.32,.36,.39,.45,.5,.56,.63,.7,.79,.88,1};
  	const int Nbins = 10;
  	TH1F *fathist = new TH1F("fathist","Xj fat based quench", Nbins, bins);
  	fathist->SetXTitle("Xj ratio (Pt2/Pt1)");
  	fathist->SetYTitle("Count");
  	dijet_tree->Draw("XjR>>fathist","","goff");
  	fathist->Sumw2(kTRUE);
  	fathist->Scale(1/fathist->Integral());
  	fathist->SetMarkerStyle(20);
  	fathist->SetMarkerColor(kRed);
  	fathist->SetLineColor(kRed);

  	TH1F *fathistc = new TH1F("fathistc","Xj fat based quench", Nbins, bins);
  	fathistc->SetXTitle("Xj ratio (Pt2/Pt1)");
  	fathistc->SetYTitle("Count");
  	dijet_tree->Draw("Xj>>fathistc","","goff");
  	fathistc->Sumw2(kTRUE);
  	fathistc->Scale(1/fathistc->Integral());
  	fathistc->SetMarkerStyle(20);

  	TLegend *tl = new TLegend(.4,.7,.6,.9);
  	tl->AddEntry(fathist, "R-quench", "l");
  	tl->AddEntry(fathistc, "control", "l");

  	fathistc->Draw("P");
  	fathist->Draw("same P");
  	tl->Draw();
  	tc->Print("fatrate 101.pdf");

  	int nRbins = 13;
  	float Rbins[] = {0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65};
  	/* make a log Z*/
  	TH2F *XR = new TH2F("XR", "Xj wrt radius",nRbins, Rbins,Nbins,bins);
  	XR->SetXTitle("radius of wider jet");
  	XR->SetYTitle("Xj ratio");
  	dijet_tree->Draw("Xj:LeadR0>>XR","","goff");
  	tc->SetLogz();
  	XR->Draw("colz");
  	std::string printfile = "XjwrtR"+outfile;
  	tc->Print(printfile.c_str());

  	TH1F *Jr = new TH1F("Jr","radius of wider jet", nRbins,Rbins);
  	Jr->SetXTitle("radius of wider jet");
  	Jr->SetYTitle("count");
  	dijet_tree->Draw("LeadR0>>Jr","","goff");
  	Jr->Draw();
  	printfile = "JetR"+outfile;
  	tc->Print(printfile.c_str());
}