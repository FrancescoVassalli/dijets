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
    gStyle->SetErrorX(0);
  	dijet_tree->Add(infile.c_str());

  	Double_t bins[] = {.32,.36,.39,.45,.5,.56,.63,.7,.79,.88,1};
  	const int Nbins = 10;

  	int nRbins = 12;
  	Double_t Rbins[] = {.05,.15,.25,.35,.45,.55,.65,.75,.85,.95,1.05,1.15};
  	TH2F *XR = new TH2F("XR", "Xj wrt radius",nRbins-1,Rbins,Nbins,bins);
  	XR->SetXTitle("radius of wider jet");
  	XR->SetYTitle("Xj ratio");
  	dijet_tree->Draw("Xj:LeadR0>>XR","","goff");
  	tc->SetLogz();
  	XR->Draw("colz");
  	std::string printfile = "XjwrtR"+outfile;
  	tc->Print(printfile.c_str());

    std::string selectname = "XjQ";
    std::string temp;

    int quenchN = 3;
    TH1F *qH[quenchN];
    for(int i=0; i<quenchN; i++){
      temp = selectname+std::to_string(i);
      qH[i]= new TH1F("qH", temp.c_str(),Nbins,bins);
      qH[i]->SetXTitle("Xj ratio");
      qH[i]->SetYTitle("count");
      temp = temp +">>qH";
      dijet_tree->Draw(temp.c_str(),"","goff");
      qH[i]->Sumw2(kTRUE);
      if(qH[i]->Integral()!=0)
        qH[i]->Scale(1/qH[i]->Integral(),"width");
      qH[i]->SetMarkerStyle(20);
      qH[i]->SetLineWidth(3);
      qH[i]->Draw("P");
      printfile = selectname+std::to_string(i)+outfile;
      tc->Print(printfile.c_str());
    }


  	TH1F *Jr = new TH1F("Jr","radius of wider jet",40,0,1.2);
  	Jr->SetXTitle("radius of wider jet");
  	Jr->SetYTitle("count");
  	dijet_tree->Draw("LeadR0>>Jr","","goff");
  	Jr->Draw();
  	printfile = "JetR"+outfile;
  	tc->Print(printfile.c_str());

    TH1F *rselect[nRbins];
    selectname = "XjR";
    float width;
    TH1F *htemp = new TH1F("temp","",Nbins,bins);
    for(int i=1; i<nRbins-1;i++){
      temp = selectname+std::to_string(i*.1);
      rselect[i]= new TH1F(Form("h%i",i),temp.c_str(),Nbins,bins);
      rselect[i]->SetXTitle("Xj ratio");
      rselect[i]->SetYTitle("count");
      temp = selectname+std::to_string(2*i)+">>h"+std::to_string(i);
      dijet_tree->Draw(temp.c_str(),"","goff");
      temp = selectname+std::to_string(2*i+1)+">>temp";
      dijet_tree->Draw(temp.c_str(),"","goff");
      rselect[i]->Sumw2(kTRUE);
      htemp->Sumw2(kTRUE);
      rselect[i]->Add(htemp,1);
      if(rselect[i]->Integral()!=0)
        rselect[i]->Scale(1/rselect[i]->Integral(),"width");
      rselect[i]->SetMarkerStyle(20);
      rselect[i]->SetLineWidth(3);
      rselect[i]->Draw("P");
      printfile = selectname+std::to_string(i)+outfile;
      tc->Print(printfile.c_str());
    }

    delete tc;
    tc = new TCanvas("tc","1D plots", 1700, 1000);
    tc->Divide(5,2,.005,.01);

    TLatex lax; 
    for(int i=1; i<nRbins-1;i++){
      tc->cd(i);
      tc->SetGrid();
      lax.SetTextSize(.07);
      temp = selectname+std::to_string(i*.1);
      lax.SetTextAlign(11);
      lax.DrawLatex(.1,.1,temp.c_str());
      //lax[i-1].PaintLatex(.3,.3,0,.5,temp.c_str());
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(.2);
      gPad->SetTicky(1);
      rselect[i]->SetTitleSize(.07);
      rselect[i]->SetYTitle("");
      rselect[i]->GetXaxis()->SetLabelSize(.06);
      rselect[i]->GetYaxis()->SetLabelSize(.06);
      rselect[i]->Draw("P");
    }
    printfile = selectname+outfile;
    tc->SaveAs(printfile.c_str());
}