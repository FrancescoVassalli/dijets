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


void FatHist(int fileN){
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
  	int filecount = 1;
  	std::string filebegin = "dijet";
  	std::string fileend = "001.root";
  	std::string infile = filebegin+std::to_string(filecount)+fileend;
  	gStyle->SetOptStat(0);
  	while(filecount<=fileN){
  		dijet_tree->Add(infile.c_str());
  		filecount++;
  		infile = filebegin+std::to_string(filecount)+fileend;
  	}

  	float bins[] = {.32,.36,.39,.45,.5,.56,.63,.7,.79,.88,1};
  	const int Nbins = 10;
  	TH1F *fathist = new TH1F("fathist","Xj fat based quench", Nbins, bins);
  	fathist->SetXTitle("Xj ratio (Pt2/Pt1)");
  	fathist->SetYTitle("Count");
  	dijet_tree->Draw("X1>>fathist","","goff");
  	fathist->Sumw2(kTRUE);
  	fathist->Scale(1/fathist->Integral());
  	fathist->SetMarkerStyle(20);
  	fathist->SetMarkerColor(kRed);
  	fathist->SetLineColor(kRed);

  	TH1F *fathist2 = new TH1F("fathist2","Xj fat based quench", Nbins, bins);
  	fathist2->SetXTitle("Xj ratio (Pt2/Pt1)");
  	fathist2->SetYTitle("Count");
  	dijet_tree->Draw("X2>>fathist2","","goff");
  	fathist2->Sumw2(kTRUE);
  	fathist2->Scale(1/fathist2->Integral());
  	fathist2->SetMarkerStyle(20);
  	fathist2->SetMarkerColor(kOrange);
  	fathist2->SetLineColor(kOrange);

  	TH1F *fathistc = new TH1F("fathistc","Xj fat based quench", Nbins, bins);
  	fathistc->SetXTitle("Xj ratio (Pt2/Pt1)");
  	fathistc->SetYTitle("Count");
  	dijet_tree->Draw("Xj>>fathistc","","goff");
  	fathistc->Sumw2(kTRUE);
  	fathistc->Scale(1/fathistc->Integral());
  	fathistc->SetMarkerStyle(20);

  	TLegend *tl = new TLegend(.4,.7,.6,.9);
  	tl->AddEntry(fathist, "antilinear", "l");
  	tl->AddEntry(fathist2, "antiquadratic","l");
  	tl->AddEntry(fathistc, "control", "l");

  	fathistc->Draw("P");
  	fathist->Draw("same P");
  	fathist2->Draw("same P");
  	tl->Draw();
  	tc->Print("fatrate 101.pdf");

  	for(int i=0; i<Nbins; i++){
  		std::cout<<fathist2->GetBinContent(i)<<'\n';
  	}
  	
}