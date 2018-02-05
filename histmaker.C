#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>
//ajust by bin width

/*
	calculates the Chi^2 value for MC compared to experimental data
*/
double myChi(TH1F *h, TH1F *h2, int NBINs){
	double sum=0;
	double data;
	double model;
	double error;
	for(int i=0; i<NBINs;i++){
		model = h->GetBinContent(i);
		data= h2->GetBinContent(i);
		error= h->GetBinError(i);
		std::cout<<"model: "<<model<<" data: "<<data<<" error: "<<error<<'\n';
		if(model!=0&&data!=0&&error!=0)
			sum+= TMath::Power((model-data),2)/TMath::Power(error,2);
	}
	return sum;
}
double myChi(TH1F *h, TH1F *h2, int NBINs, int nFree){
	double sum=0;
	double data;
	double model;
	double error;
	for(int i=0; i<NBINs;i++){
		model = h2->GetBinContent(i);
		data= h->GetBinContent(i);
		error= h2->GetBinError(i);
		if(model!=0&&data!=0&&error!=0)
			sum+= TMath::Power((model-data),2)/TMath::Power(error,2);
	}
	return sum/nFree;
}

void histmaker(int fileN)
{
  TCanvas *tc = new TCanvas();
  bool lowpT = true;
  const int DFreedom =10;
  std::string filebegin = "dijet";
  std::string fileend = "001.root";
  float stateMax;
  std::string file;
  if(lowpT){
    stateMax=3.5;
    file = " 100.pdf";
  }
  else{
    stateMax = 4.5;
    file = " 200.pdf";
  }
  TChain *dijet_tree = new TChain("tree100");
  int filecount = 1;
  std::string filename = filebegin+std::to_string(filecount)+fileend;
  while(filecount<=fileN){
    dijet_tree->Add(filename.c_str());
    filecount++;
    filename = filebegin+std::to_string(fileN)+fileend;
  }
  TLegend *tl=new TLegend(.4,.7,.6,.9);
  const int HISTN = 55; ///         0   1      2       3    4   5    6       7    8      9    10    11    12    13    14    15    16  17      18    19     20    21  22  23    24  25    26   27     28    29   30      31     32    33       34    35     36    37     38      39    40     41       42     43     44          45     46    47   48   49  50   51    52        53		54
  std::string histnames[HISTN] = {"Xj","QQ1","QQ2","QQ3","QQA","QQB","QG1","QG2","QG3","QGA","QGB","GQ1","GQ2","GQ3","GQA","GQB","GG1","GG2","GG3","GGA","GGB","X1","X2","X3","XA","XB","XC","QQC", "QGC","GQC","GGC","RAQQ","RAQG","RAGQ","RAGG","RBQQ","RBQG","RBGQ","RBGG","RCQQ","RCQG","RCGQ","RCGG","xrate","quadFat", "linFat","ZR3","Z4","XD","XE","X4","X5","xrateB","xrateC", "XP"};
  std::string treeRoute;
  std::string Qtemp;
  TH1F *h[HISTN] = {NULL};
  TH1F *hq[HISTN]={NULL};
  gStyle->SetOptStat(0);
  std::string counts[HISTN];
  std::string temp;
  float width;
  float bins[] = {.32,.36,.39,.45,.5,.56,.63,.7,.79,.88,1};
  const int BINN =10;
  float atlasV[]={.05,.3,.8,1.8,2.55,2.35,1.78,1.4,1.35,1.33};
  float atlasE[]={.05,.2,.35,.42,.45,.25,.2,.2,.2,.2};
  float atlasV2[]={.2,.45,.6,.8,1,1.2,1.5,1.8,2.22,2.5};
  float atlasE2[]={.05,.1,.12,.15,.1,.1,.1,.1,.1};
  TH1F *hAtlas = new TH1F("Atlas","Atlas",10,bins); // the 100 pT
  TH1F *hAtlas2 = new TH1F("Atlas2","Atlas2",10,bins); // outer centrality
  TH1F *QAtlas = new TH1F("QAtlas","QAtlas",10,bins);
  hAtlas->SetXTitle("Xj ratio (Pt2/Pt1)");
  hAtlas->SetYTitle("Count");
  hAtlas2->SetXTitle("Xj ratio (Pt2/Pt1)");
  hAtlas2->SetYTitle("Count");
  	for(int i=0; i<BINN;i++){
   	 	hAtlas->SetBinContent(i+1,atlasV[i]);
    	hAtlas->SetBinError(i+1,atlasE[i]);
    	hAtlas2->SetBinContent(i+1,atlasV2[i]);
    	hAtlas2->SetBinError(i+1,atlasE2[i]);
  	}
  QAtlas->Add(hAtlas,hAtlas2,1,-1);
  hAtlas->Draw();
  hAtlas2->SetLineColor(kRed);
  QAtlas->SetLineColor(kGreen+2);
  tl->AddEntry(hAtlas,"Atlas 100","l");
  hAtlas2->Draw("same");
  QAtlas->Draw("same");
  tl->AddEntry(hAtlas2,"Atlas 200","l");
  tl->AddEntry(QAtlas, "Quench","l");
  tl->Draw("same");
  tc->Print("AtlasRoot.pdf");
  tc->Clear("D");
  tl->Clear("D");
  for(int i=0; i<HISTN; i++){
    h[i] = new TH1F(histnames[i].c_str(),histnames[i].c_str(),10, bins);
    Qtemp = histnames[i]+"quench";
    hq[i] = new TH1F(Qtemp.c_str(),Qtemp.c_str(),10,bins);
    h[i]->SetXTitle("Xj ratio (Pt2/Pt1)");
    h[i]->SetYTitle("Count");
    treeRoute = histnames[i]+">>"+histnames[i];
    dijet_tree->Draw(treeRoute.c_str(),"","goff");
    counts[i] = std::to_string(h[i]->Integral());
    h[i]->Scale(1/h[i]->Integral());
    for(int j=0; j<BINN;j++){
      width = bins[j+1]-bins[j];
      h[i]->SetBinContent(j+1,h[i]->GetBinContent(j+1)/width);
    }
    h[i]->SetMarkerStyle(20);
    h[i]->SetLineWidth(3);
  }
  h[0]->Draw("P");
  hAtlas->Draw("same P");
  tl->AddEntry(h[0], "Pythia Control","l");
  tl->AddEntry(hAtlas, "Atlas","l");
  temp = "Xj control"+file;
  tc->Print(temp.c_str());
  tc->Clear("D");
  tl->Clear("D");


  hAtlas->Draw("P");
  h[21]->Draw("same P");
  h[22]->Draw("same P");
  h[23]->Draw("same P");
  h[50]->Draw("same P");
  h[51]->Draw("same P");
  hAtlas->SetMaximum(3);
  hAtlas->SetMinimum(0);
  h[22]->SetLineColor(kGreen+2);
  h[23]->SetLineColor(kRed);
  h[50]->SetLineColor(kOrange);
  h[51]->SetLineColor(41);

  delete tl;
  tl=new TLegend(.4,.1,.7,.35);

  Double_t res[BINN];
  double compare = myChi(h[21], hAtlas,BINN, DFreedom);
  double chiT = h[21]->Chi2Test(hAtlas,"WW CHI2", res)/DFreedom;
  temp = "Ab(0,20) Chi2: "+std::to_string(chiT)+" - "+std::to_string(compare);
  tl->AddEntry(h[21],temp.c_str(),"l");

  compare = myChi(h[22], hAtlas, BINN, DFreedom);
  chiT = h[22]->Chi2Test(hAtlas,"WW CHI2", res)/DFreedom;
  temp = "Ab(10,20) Chi2: "+std::to_string(chiT)+" - "+std::to_string(compare);
  tl->AddEntry(h[22],temp.c_str(),"l");

  compare = myChi(h[23], hAtlas, BINN, DFreedom);
  chiT = h[23]->Chi2Test(hAtlas,"WW CHI2", res)/DFreedom;
  temp = "Ab(20,20) Chi2: "+std::to_string(chiT)+" - "+std::to_string(compare);
  tl->AddEntry(h[23],temp.c_str(),"l");

  compare = myChi(h[50], hAtlas, BINN, DFreedom);
  chiT = h[50]->Chi2Test(hAtlas,"WW CHI2", res)/DFreedom;
  temp = "Ab(20,30) Chi2: "+std::to_string(chiT)+" - "+std::to_string(compare);
  tl->AddEntry(h[50],temp.c_str(),"l");

  compare = myChi(h[51], hAtlas, BINN, DFreedom);
  chiT = h[51]->Chi2Test(hAtlas,"WW CHI2", res)/DFreedom;
  temp = "Ab(20,40) Chi2: "+std::to_string(chiT)+" - "+std::to_string(compare);
  tl->AddEntry(h[51],temp.c_str(),"l");
  tl->AddEntry(hAtlas, "Atlas","l");
  tl->Draw();
  tc->Print("AB fit.pdf");
  tc->Clear("D");
  tl->Clear();

  hAtlas->Draw("P");
  h[24]->Draw("same P");
  h[25]->Draw("same hist");
  h[26]->Draw("same hist");
  h[48]->Draw("same hist");
  h[49]->Draw("same hist");
  hAtlas->SetMaximum(3);
  hAtlas->SetMinimum(0);
  h[25]->SetLineColor(kGreen+2);
  h[26]->SetLineColor(kRed);
  h[48]->SetLineColor(kOrange);
  h[49]->SetLineColor(41);
  chiT=hAtlas->Chi2Test(h[24],"WW CHI2",res)/DFreedom;
  temp="Rel(.2,.05) Chi2: "+std::to_string(chiT);
  tl->AddEntry(h[24],temp.c_str(),"l");
  chiT=hAtlas->Chi2Test(h[25],"WW CHI2",res)/DFreedom;
  temp="Rel(.2,.1) Chi2: "+std::to_string(chiT);
  tl->AddEntry(h[25],temp.c_str(),"l");
  chiT=hAtlas->Chi2Test(h[26],"WW CHI2",res)/DFreedom;
  temp="Rel(.2,.2)" +std::to_string(chiT);
  tl->AddEntry(h[26],temp.c_str(),"l");
  chiT=hAtlas->Chi2Test(h[48],"WW CHI2",res)/DFreedom;
  temp="Rel(.3,.2)" +std::to_string(chiT);
  chiT=hAtlas->Chi2Test(h[49],"WW CHI2",res)/DFreedom;
  temp="Rel(.1,.2)" +std::to_string(chiT);
  tl->AddEntry(h[48],temp.c_str(),"l");
  tl->AddEntry(hAtlas, "Atlas","l");
  tl->Draw();
  tc->Print("Rel fit.pdf");
  tc->Clear("D");
  tl->Clear("D");

  h[0]->SetMaximum(5);
  h[0]->Draw();
  h[1]->Draw("same hist");
  h[6]->Draw("same hist");
  h[11]->Draw("same hist");
  h[16]->Draw("same hist");
  h[1]->SetLineColor(kRed);
  h[6]->SetLineColor(kOrange);
  h[11]->SetLineColor(41);
  h[16]->SetLineColor(kCyan);
  delete tl;
  tl=new TLegend(.4,.7,.6,.9);
  tl->AddEntry(h[0],"control","l");
  tl->AddEntry(h[1],"QQ unquenched","l");
  tl->AddEntry(h[6],"QG","l");
  tl->AddEntry(h[11],"GQ","l");
  tl->AddEntry(h[16],"GG","l");
  tl->Draw();
  temp = "unquenched flavor" + file;
  tc->Print(temp.c_str());
  tc->Clear("D");
  tl->Clear("D");

  h[0]->SetMaximum(stateMax);
  h[0]->Draw();
  h[23]->Draw("same hist");
  h[50]->Draw("same hist");
  h[51]->Draw("same hist");
  tl->AddEntry(h[0],"control","l");
  tl->AddEntry(h[23],"GAUS(20,20)","l");
  tl->AddEntry(h[50],"GAUS(20,30)","l");
  tl->AddEntry(h[51],"GAUS(20,40)","l");
  tl->Draw();
  temp = "Abs large" + file;
  tc->Print(temp.c_str());
  tc->Clear("D");
  tl->Clear("D");

  h[0]->Draw();
  tl->AddEntry(h[0],"control","l");
  h[24]->Draw("same");
  h[24]->SetLineColor(41);
  h[25]->Draw("same hist");
  h[26]->Draw("same hist");
  h[49]->Draw("same hist");
  tl->AddEntry(h[24],"GAUS(.2,.05)","l");
  tl->AddEntry(h[25],"GAUS(.2,.1)","l");
  tl->AddEntry(h[26],"GAUS(.2,.2)","l");
  tl->AddEntry(h[49],"GAUS(.1,.2)","l");
  tl->Draw();
  temp = "Relative"+file;
  tc->Print(temp.c_str());
  tl->Clear("D");
  tc->Clear("D");

  h[0]->Draw();
  h[43]->SetLineColor(kGreen+2);
  h[52]->SetLineColor(kRed);
  h[53]->SetLineColor(41);
  tl->AddEntry(h[0],"control","l");
  h[43]->Draw("same hist");
  h[52]->Draw("same hist");
  h[53]->Draw("same hist");
  tl->AddEntry(h[43],"1.5 GeV per Particle","l");
  tl->AddEntry(h[52],"2 per Particle","l");
  tl->AddEntry(h[53],"1 per Particle","l");
  tl->Draw();
  temp = "Multiplicity"+file;
  tc->Print(temp.c_str());
  tl->Clear("D");
  tc->Clear("D");

  h[23]->SetMaximum(stateMax);
  h[23]->Draw();
  temp  = "GAUS(20,20) " + counts[23];
  tl->AddEntry(h[23],temp.c_str(),"l");
  h[23]->SetLineColor(kGreen+2);
  h[3]->SetLineColor(kRed);
  h[8]->SetLineColor(kOrange);
  h[18]->SetLineColor(41);
  h[13]->SetLineColor(kCyan);
  h[3]->Draw("same hist");
  temp = "QQ " + counts[3];
  tl->AddEntry(h[3],temp.c_str(),"l");
  h[8]->Draw("same");
  temp = "QG " + counts[8];
  tl->AddEntry(h[8],temp.c_str(),"l");
  h[13]->Draw("same hist");
  temp = "GQ " + counts[13];
  tl->AddEntry(h[13],temp.c_str(),"l");
  h[18]->Draw("same hist");
  temp = "GG " + counts[18];
  tl->AddEntry(h[18],temp.c_str(),"l");
  temp = "Abs Parton"+file;
  tl->Draw();
  tc->Print(temp.c_str());
  tl->Clear("D");
  tc->Clear("D");

  h[26]->SetMaximum(stateMax);
  h[26]->Draw();
  tl->AddEntry(h[26],"GAUS(.2,.2)","l");
  h[26]->SetLineColor(kGreen+2);
  h[27]->SetLineColor(kRed);
  h[28]->SetLineColor(kOrange);
  h[29]->SetLineColor(41);
  h[30]->SetLineColor(kCyan);
  h[27]->Draw("same hist");
  tl->AddEntry(h[27],"QQ","l");
  h[28]->Draw("same hist");
  tl->AddEntry(h[28],"QG","l");
  h[29]->Draw("same hist");
  tl->AddEntry(h[29],"GQ","l");
  h[30]->Draw("same hist");
  tl->AddEntry(h[30],"GG","l");
  temp = "Rel Parton"+file;
  tl->Draw();
  tc->Print(temp.c_str());
  tl->Clear("D");
  tc->Clear("D");

  h[43]->SetMaximum(stateMax);
  h[43]->Draw();
  tl->AddEntry(h[43],"1.5 GeV per Particle","l");
  h[43]->SetLineColor(kGreen+2);
  h[31]->SetLineColor(kRed);
  h[32]->SetLineColor(kOrange);
  h[33]->SetLineColor(41);
  h[34]->SetLineColor(kCyan);
  h[31]->Draw("same hist");
  tl->AddEntry(h[31],"QQ","l");
  h[32]->Draw("same hist");
  tl->AddEntry(h[32],"QG","l");
  h[33]->Draw("same hist");
  tl->AddEntry(h[33],"GQ","l");
  h[34]->Draw("same hist");
  tl->AddEntry(h[34],"GG","l");
  temp = "Mult Parton"+file;
  tl->Draw();
  tc->Print(temp.c_str());
  tl->Clear("D");

  h[54]->SetMaximum(stateMax);
  h[54]->Draw();
  tl->AddEntry(h[54],"Inerpolation", "l");
  temp = "Interpolation"+file;
  tc->Print(temp.c_str());
  tl->Clear("D");

  tc->SetLogy();
  TH1F *h_jet = new TH1F("jets", "jets", 88, 0,450);
  dijet_tree->Draw("highjet>>jets","","goff");
  h_jet->Draw();
  tc->Update();
  TLine *myLine = new TLine(100, 1, 100, 100000); 
  myLine->SetLineColor(kRed);
  myLine->Draw();
  temp = "highjets"+file;
  tc->Print(temp.c_str());
  tc->Clear("D");

  TH1F *ljet = new TH1F("ljets", "ljets", 44, 0,450);
  dijet_tree->Draw("lowjet>>ljets","","goff");
  ljet->Draw();
  temp = "lowjets"+file;
  tc->Print(temp.c_str());
}
