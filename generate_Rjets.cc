/*
parton pairing isnt working because the eType function times out
*/

#include "Pythia8/Pythia.h"
using namespace Pythia8;
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "spline.h"
#include <fstream>
#include <array>

Double_t E= 2.71828182845904523536;

struct Jet{
  float pT, phi, y ,r, fatratio;
  int mult;
};
struct Parton{
  float px, py, pz, eta;
  int id;
};
struct XjT{
  float Xj;
  float fat;
  short type;
};

inline float randomPositive(double mean, double sig){
  float r = gRandom->Gaus(mean, sig);
  while(r<0){
    r = gRandom->Gaus(mean, sig);
  }
  return r;
}

float fixXj(float *j1, float *j2, int min, int max){
  if(*j1<*j2){
    float temp = *j1;
    *j1 = *j2;
    *j2 = temp;
  }
  if(*j1<=0||*j2<=0) 
    return 0;
  if(*j1<min||*j2<25||*j1>max)
    return 0;
  return *j2/(*j1);
}
float fixXj(float *j1, float *j2, int min){
  if(*j1<*j2){
    float temp = *j1;
    *j1 = *j2;
    *j2 = temp;
  }
  if(*j1<=0||*j2<=0)
    return 0;
  if(*j1<min||*j2<25)
    return 0;
  return *j2/(*j1);
}
inline float fixXj(float j1, float j2, int min){
  if(j1<min||j2<25)
    return 0;
  return j2/(j1);
}
inline float fixXj(float j1, float j2, int min, int max){
  if(j1<min||j2<25||j1>max||j2>max)
    return 0;
  return j2/(j1);
}

float delPhi(float phi1, Parton p1){
  float r = TMath::Abs(phi1 - TMath::ATan2(p1.py,p1.px));
  while(r>TMath::Pi())
    r-=TMath::Pi();
  return r;
}
float delphi(float phi1, float phi2){
  float r = TMath::Abs(phi1-phi2);
  while(r>TMath::Pi())
    r-=TMath::Pi();
  return r;
}

inline float quadrature(float x, float y , float z){
  return TMath::Sqrt(x*x+y*y+z*z);
}


float etaToRad(float eta1){
  Double_t e = (double) eta1*-1;
  float r= 2* TMath::ATan2((double)TMath::Power(E,e),1);
  while(r>2*TMath::Pi())
    r-=2*TMath::Pi();
  return r;
}

float delR(float phi1, float y1, Parton p1){
  y1 = etaToRad(y1);
  float y2 = etaToRad(p1.eta);
  phi1 = delPhi(phi1, p1);
  return quadrature(TMath::Abs(y1-y2),phi1,0);
}

inline float delR(Jet skinny, Jet fatty){
  skinny.y = etaToRad(skinny.y);
  fatty.y = etaToRad(fatty.y);
  return quadrature(TMath::Abs(skinny.y-fatty.y), delphi(skinny.phi,fatty.phi),0);
}

inline void swapPart(Parton *p1, Parton *p2){
  Parton temp = *p1;
  *p1 = *p2;
  *p2=temp;
}

inline void setLead(Parton *p1, Parton *p2, float phi, float y){
  if(delR(phi,y,*p1)>delR(phi,y,*p2))
    {
      swapPart(p1, p2);
    }
}
inline bool isQuark(Parton p1){
  if(TMath::Abs(p1.id)>=1&& TMath::Abs(p1.id)<=8)
    return true;
  return false;
}

short eType(Parton *p1,Parton *p2, float phi, float y){
  setLead(p1,p2,phi,y);
  if(isQuark(*p1)&&isQuark(*p2))return 1;
  if(p2->id==21&&p1->id==21)
    return 4;

  if( p1->id==21&&isQuark(*p2))
    {
	return 3;  //leading jet is a gluon and sub is quark
    }
  else{
    if(isQuark(*p1)&&p2->id==21) 
      return 2; //if lead is quark and sub is gluon
  }
  return -1;
}


int jetMax1(std::vector<Jet> jets, float* jet, float max){
  *jet=0;
  int r=-1;
  if(max==0){
      for(unsigned i=0; i<jets.size();i++){
	       if(jets[i].pT>*jet){
	       *jet=jets[i].pT;
         r=i;
	       }
      }
  }
  else{
    for(unsigned i=0; i<jets.size();i++){
	     if(jets[i].pT >*jet && jets[i].pT<max){
	       *jet=jets[i].pT;
         r=i;
	     }
    }
  }
  return r;
}

void jetMax2(int count, int *index, float *jet, float jet_a[], float jet_max, int min){
  *jet=0;
  *index =-1;
  if(jet_max==0){
      for(int i=0; i<count;i++){
	if(jet_a[i]>*jet){
	  *jet=jet_a[i];
	  *index=i;
	}
      }
  }
  else{
    for(int i=0; i<count;i++){
	if(jet_a[i] >*jet && jet_a[i]<jet_max && jet_a[i]>min){
	  *jet=jet_a[i];
	  *index=i;
	}
    }
  }
}

int minimum(std::vector<float> v){
  float t= v[0];
  int r=0;
  for(unsigned i=1; i<v.size(); i++){
    if(v[i]<t){
      t= v[i];
      r=i;
    }
  }
  return r;
}

 std::vector<Jet> makeFatRatio (std::vector<Jet> jet, std::vector<Jet> fat){
  std::vector<float> dR(fat.size());
  std::vector<Jet> out(jet.size());
  int mindex;
  for(unsigned i=0; i<jet.size();i++){
    for(unsigned j=0; j<fat.size();j++){
      dR[j] = delR(jet[i], fat[j]);
    }
    mindex = minimum(dR);
    out[i].fatratio = jet[i].pT/fat[mindex].pT;
 	out[i].pT =jet[i].pT;
 	out[i].y = jet[i].y;
 	out[i].phi = jet[i].phi;
 	out[i].mult = jet[i].mult;
  }
  return out;
}

void setJetData(int index, float* phi, float* y, float jet_phi[], float jet_y[]){
  if(index==-1){
    *phi=0;
    *y=0;
  }
  *phi = jet_phi[index];
  *y = jet_y[index];
}

std::vector<Jet> findR(std::vector<Jet> *jets, int nR){
	std::vector<Jet> out(nR);
	std::vector<Jet> tempout;
	for(int i=0; i<nR-1; i++){
		tempout= makeFatRatio(jets[i],jets[i+1]);
		for(unsigned int j=0; j<nR;j++){
			if(tempout[j].fatratio>.9||i==nR-2){
				out[j]=tempout[j];
				out[j].r = .9-i*.05;
			}
		}
	}
	return out;
}

void makedata(std::string filename,int fitNUM, int fitMAX, bool lowpT, int nEvent){
  TFile* f = new TFile(filename.c_str(),"RECREATE");
  TTree* t=new TTree("tree100","100pThat events");
  t->SetAutoFlush(-70000);
  std::vector<double> interE = {50,70,90,100,110,126,140,170,200}; 
  std::vector<double> interM={.5,.43,.35,.3,.25,.2,.15,.1,.05};

  Pythia pythia;
  pythia.readString("Beams:eCM = 2760.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random::setSeed = on");
  pythia.readString("Random::seed =0");
  if(lowpT)
    pythia.readString("PhaseSpace:pTHatMin = 80.");
  else
    pythia.readString("PhaseSpace:pTHatMin = 150.");
  pythia.init();
  const int nR=7;
  std::vector<SlowJet*> kT(nR);
  SlowJet *tempjet;

  std::vector<float> fatratio;
  Parton p1; Parton p2;

  const int nXj = 2;
  XjT Xjs[nXj];
  float QQ1,QG1,GQ1,GG1;
  t->Branch("Xj", &Xjs[0].Xj);
  t->Branch("XjR", &Xjs[1].Xj);

  const int nfjets=nXj*2;
  float fjets[nfjets];
  int eventType[nXj];
  int iAlag=0;
  int ipt1=0;
  std::vector<Jet> jets[nR];
  std::vector<Jet> myJets(nR);
  std::vector<Jet> tempjets(nR);
  //std::vector<Jet> tempjets[nR];
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next())
      continue;
    //std::cout<<"event loop"<<iEvent<<'\n';
    p1.px = pythia.event[5].px();
    p1.py=pythia.event[5].py();
    p1.eta = pythia.event[5].eta();
    p2.eta = pythia.event[6].eta();
    p2.px=pythia.event[6].px();
    p2.py=pythia.event[6].py();
    p1.id= pythia.event[5].id();
    p2.id= pythia.event[6].id();
    for(int i=0; i<nR; i++){  //make an array of the jets
    	tempjet= new SlowJet(-1,0.9-.05*i, 10,4,2,1);
    	tempjet->analyze(pythia.event);
    	jets[i].resize(tempjet->sizeJet());
    	for(int j=0; j<tempjet->sizeJet();j++){
    		jets[i][j].pT= tempjet->pT(i);
    		jets[i][j].y= tempjet->y(i);
    		jets[i][j].phi= tempjet->phi(i);
    		jets[i][j].mult= tempjet->multiplicity(i);
    	}
    	delete tempjet;
    }

    myJets = findR(jets, nR);// turn the 2D array into a 1D array of 90% jets
   	tempjets=myJets;
    QQ1=0;QG1=0;GQ1=0;GG1=0;
    iAlag=0;
    ipt1=0;
    jetMax1(myJets,&fjets[0],0); //1
    jetMax1(myJets,&fjets[1],fjets[0]);
//2
   /* for(int i=0; i<nR;++i){
       myJets[i].pT=myJets[i].pT-randomPositive(20,10)/(1-myJets[i].r);
     }*/
       ipt1 = jetMax1(myJets, &fjets[2],0);
       jetMax1(myJets, &fjets[3],fjets[2]);
       //eventType[iAlag++] = eType(&p1,&p2,myJets[ipt1].phi,myJets[ipt1].y);
       myJets=tempjets;
//3
    for(int i=0; i<nXj; i++){
      Xjs[i].Xj = fixXj(&fjets[2*i+1],&fjets[2*i],fitNUM,fitMAX);
      Xjs[i].type = eventType[i];
      if(Xjs[i].Xj>1){
         Xjs[i].Xj = 1/Xjs[i].Xj;
         if(Xjs[i].type ==2)
             Xjs[i].type=3;
         else if(Xjs[i].type==3)
             Xjs[i].type=2;
      }
    }
    switch(Xjs[1].type){
    case 1:
       QQ1 = Xjs[1].Xj;
       break;
    case 2:
       QG1 = Xjs[1].Xj;
       break;
    case 3:
       GQ1 = Xjs[1].Xj;
       break;
    case 4:
       GG1 = Xjs[1].Xj;
       break;
    }

    t->Fill();
  }
  t->Write();
  f->Write();
  f->Close();
  delete f;
}

int main(int argc, char *argv[]){
	std::string filename;
	int fitNUM, fitMAX;
	bool lowpT;
	int nEvent = 10000;
	if(argc!=3){
		std::cout<<"accepts 2 arguments: 1. outfile 2. low or high pT"<<'\n';
		return 1;
	}
	else{
		filename=argv[1];
		std::string temp(argv[2]);
		if(temp=="low"||temp=="l")
			lowpT=true;
		else if(temp=="high"||temp=="h")
			lowpT=false;
		else{
			std::cout<<"Incorrect pT setting"<<'\n';
			return 2;
		}
		if(lowpT){
			fitNUM =100;
    		fitMAX = 126;
		}
		else{
			fitNUM=200;
    		fitMAX=3000;
		}
		makedata(filename, fitNUM, fitMAX,lowpT, nEvent);
	}
	return 0;
}