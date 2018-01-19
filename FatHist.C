#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "TChain.h"

/**
add 1D histrograms seperated by R

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
  		infile= "Rjet21";
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

  	int nRbins = 20;
  	float Rbins[] = {0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1};
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
/*
    TH1F *rselect;
    std::string selectname = "XjR";
    std::string temp;
    for(int i=0; i<nRbins+1;i++){
      temp = selectname+std::to_string(i*.05+.05);
      rselect= new TH1F("rselect",temp.c_str(),Nbins,bins);
      rselect->SetXTitle("Xj ratio");
      rselect->SetYTitle("count");
      temp = selectname+std::to_string(i)+">>rselect";
      dijet_tree->Draw(temp.c_str(),"","goff");
      rselect->Sumw2(kTRUE);
      rselect->Scale(1/rselect->Integral());
      rselect->Draw("P");
      printfile = selectname+std::to_string(i)+outfile;
      tc->Print(printfile.c_str());
      delete rselect;
    }*/
}