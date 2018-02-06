
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
#include <stdexcept>

Double_t E= 2.71828182845904523536;

struct Jet{
  float pT, phi, y;
  float r=-1;
  float fatratio=0;
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
  float r;
};

float maxFloat(float, float);
void printJets(std::vector<Jet>);
void printJets(std::vector<Jet>,std::string);
bool nullJets(std::vector<Jet> v);

float randomPositive(double mean, double sig){
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
 float fixXj(float j1, float j2, int min){
  if(j1<min||j2<25)
    return 0;
  return j2/(j1);
}
 float fixXj(float j1, float j2, int min, int max){
  if(j1<min||j2<25||j1>max||j2>max)
    return 0;
  return j2/(j1);
}
 float fixXj(float j1, float j2){
	if(j1>j2)
		return j2/j1;
	else
		return j1/j2;
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

float delR(Jet skinny, Jet fatty){
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
  if(v.size()==0){
    cout<<"invalid size in function: int minimum(std::vector<float> v)"<<std::endl;
    throw std::invalid_argument("invalid array length");
  }
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

/*
jet->size() msut be greater than 0
*/
std::vector<Jet> makeFatRatio (std::vector<Jet> jet, std::vector<Jet> fat){ // what if jet->size is zero
  std::vector<Jet> out(jet.size());
  std::vector<float> dR(fat.size());
  int mindex =-1;
  for(unsigned i=0; i<jet.size();i++){  
    for(unsigned j=0; j<fat.size();j++){  // find the differnce in radius between the current jet and all the larger jets
      dR[j] = delR(jet[i],fat[j]);
    }
    if(fat.size()>0){  //if there are bigger jets 
      mindex = minimum(dR); // locate the large jet in the same direction as the current jet
      jet[i].fatratio = jet[i].pT/fat[mindex].pT; // set their fat ratio
    }
    else{
      out[i].fatratio = 1;
      std::cout<<"Fat Size Warning"<<std::endl;
    }
 	  out[i].pT =jet[i].pT; //load the current jet into the out table
 	  out[i].y = jet[i].y;
 	  out[i].phi = jet[i].phi;
 	  out[i].mult = jet[i].mult;
    out[i].fatratio = jet[i].fatratio;
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

float findRatio(std::vector<Jet> jets, std::vector<Jet> fats){
  float minR=-1;
  std::vector<Jet> tempout;
  tempout = makeFatRatio(jets,fats);
  std::cout<<tempout.size()<<std::endl;
  for(unsigned i=0; i<tempout.size(); i++){
    if(minR<0||tempout[i].fatratio<minR)
      minR = tempout[i].fatratio;
  }
  /*
	for(int i=0; i<nR-1; i++){
		tempout= makeFatRatio(&jets[i],&jets[i+1]);
		for(unsigned j=0; j<tempout.size();j++){
			if((tempout[j].fatratio>.9||i==nR-2)&&tempout[j].pT>0){
				tempout[j].r=i*step+minR;
				//cout<<tempout[j].pT<<'\n';
				out.push_back(tempout[j]);
			}
		}
	}*/
	return minR;
}

std::vector<Jet> fillJets(std::vector<Jet> myJets, std::vector<Jet> Tjets){
  std::vector<Jet> out;
  if(!nullJets(myJets)){
    std::vector<float> dR(Tjets.size());
    int mindex;
    for(unsigned i=0; i<myJets.size(); i++){
     if(myJets[i].fatratio<.9){
      // std::cout<<"delR: ";
       for(unsigned j=0; j<Tjets.size();j++){
         dR[j] = delR(myJets[i],Tjets[j]);
        // std::cout<<dR[j]<<" ";
         //std::cout<<"loop if for"<<std::endl;
       }
       mindex = minimum(dR);/*
       std::cout<<"min: "<<mindex<<std::endl;
       out[i].pT = Tjets[mindex].pT;
       out[i].y = Tjets[mindex].y;
       out[i].phi =Tjets[mindex].phi;
       out[i].mult = Tjets[mindex].mult;*/
       out.push_back(Tjets[mindex]);
      //std::cout<<"\n";
      }
      else{
        out.push_back(myJets[i]);
      }
    }
  }
  else{
    out=Tjets;
  }
  return out;
}

void jetFilter(std::vector<Jet> *jets, int fitNUM, int fitMAX){
	for(unsigned i=0;i<jets->size();i++){
		if((*jets)[i].pT<fitNUM||(*jets)[i].pT>fitMAX)
			jets->erase(jets->begin()+i);
	}
}
inline bool nullJets(std::vector<Jet> v){
  bool r = true;
  for(unsigned i=0; i<v.size();i++){
    if(v[i].pT!=0)
      r=false;
  }
  return r;
}

inline float getMaxR(Jet j1, Jet j2){
	if(j2.r>j1.r)
		return j2.r;
	else
		return j1.r;
}
inline float maxFloat(float j1, float j2){
	if(j2>j1)
		return j2;
	else
		return j1;
}

inline void printJets(std::vector<Jet> v){
  std::cout<<"pTs: ";
  for(unsigned i=0; i<v.size();i++){
    std::cout<<v[i].pT<<" ";
  }
  std::cout<<std::endl;
}
inline void printJets(std::vector<Jet> v, std::string title){
  std::cout<<title<<" pTs: ";
  for(unsigned i=0; i<v.size();i++){
    std::cout<<v[i].pT<<" ";
  }
  std::cout<<std::endl;
}

void makedata(std::string filename,int fitNUM, int fitMAX, bool lowpT, int nEvent){
  TFile* f = new TFile(filename.c_str(),"RECREATE");
  TTree* t=new TTree("tree100","events");
  //t->SetAutoFlush(-70000);

  Pythia pythia;
  pythia.readString("Beams:eCM = 2760.");
  pythia.readString("HardQCD:all = on");
  /*
  pythia flair promptphoton:all
  for gamma jet pair
  limit pThat around 40
  rebin for gamma
  */
  pythia.readString("Random::setSeed = on");
  pythia.readString("Random::seed =0");
  if(lowpT)
    pythia.readString("PhaseSpace:pTHatMin = 80.");
  else
    pythia.readString("PhaseSpace:pTHatMin = 150.");
  pythia.init();

  const float step =.05;
  const float minR = .05;
  const int nR=24;
  std::vector<SlowJet*> tempjet(2);
  int tempjetcounter;
  
  const short preN =3;
  const int nXj = preN+nR;
  XjT Xjs[nXj];
  t->Branch("Xj", &Xjs[0].Xj); 
  t->Branch("XjQ0", &Xjs[1].Xj);
  t->Branch("XjQ1", &Xjs[2].Xj);
  t->Branch("LeadR0",&Xjs[0].r);
  std::string branchname = "XjR";
  std::string stringtemp;
  /*for(int i=0; i<nR;i++){
  	stringtemp = branchname+std::to_string(i);
  	t->Branch(stringtemp.c_str(),&Xjs[preN+i].Xj);
  }
*/
  float eventRadius;
  float eventRatio;
  const int nfjets=nXj*2;
  float fjets[nfjets];
  int ipt1=0;
  float fR[nfjets];
  std::vector<Jet> jets;
  std::vector<Jet> myJets;
  std::vector<Jet> fats;
  std::vector<Jet> tempjets;
  //Parton p1; Parton p2;
  int jetfindcounter=0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
  	if(iEvent%30==0)  
  		cout<<"Event N: "<<iEvent<<'\n';
    if (!pythia.next()){
      cout<<"pythia.next() failed"<<"\n";
      continue;
    }
    //std::cout<<"event loop"<<iEvent<<'\n';
   /* p1.px = pythia.event[5].px(); //old code for parton typing
    p1.py=pythia.event[5].py();
    p1.eta = pythia.event[5].eta();
    p2.eta = pythia.event[6].eta();
    p2.px=pythia.event[6].px();
    p2.py=pythia.event[6].py();
    p1.id= pythia.event[5].id();
    p2.id= pythia.event[6].id();*/
    /*
    
    for(int i=0; i<nR; i++){  //old code to make an array of the jets
    	tempjet= new SlowJet(-1,step*i+minR, 10,4,2,1); 
    	tempjet->analyze(pythia.event);
    	jets[i].resize(tempjet->sizeJet());
    	for(int j=0; j<tempjet->sizeJet();j++){
    		jets[i][j].pT= tempjet->pT(i);
    		//cout<<tempjet->pT(i)<<'\n';
    		jets[i][j].y= tempjet->y(i);
    		jets[i][j].phi= tempjet->phi(i);
    		jets[i][j].mult= tempjet->multiplicity(i);
    	}
    	delete tempjet;
    }*/
    tempjet[0] = new SlowJet(-1,1,10,4,2,1);  //set up the comparison jets
    tempjet[0]->analyze(pythia.event);
    fats.resize(tempjet[0]->sizeJet());
    for(unsigned i=0; i<fats.size(); i++){
      fats[i].pT = tempjet[0]->pT(i);
      fats[i].y= tempjet[0]->y(i);
      fats[i].phi= tempjet[0]->phi(i);
      fats[i].mult= tempjet[0]->multiplicity(i);
      fats[i].r = 1;
    }

    tempjet[1]= new SlowJet(-1,minR,10,4,2,1); // set up the smallest jets 
    jetfindcounter=1;
    tempjet[1]->analyze(pythia.event);
    myJets.resize(tempjet[1]->sizeJet());
    for(unsigned i=0; i<jets.size(); i++){
      myJets[i].pT = tempjet[1]->pT(i);
      myJets[i].y= tempjet[1]->y(i);
      myJets[i].phi= tempjet[1]->phi(i);
      myJets[i].mult= tempjet[1]->multiplicity(i);
      myJets[i].r = minR;
    }

    eventRatio=0;
    tempjetcounter=1;
    while(eventRatio<.9&&jetfindcounter<=19){
      tempjetcounter++;
      tempjet.push_back( new SlowJet(-1,step*jetfindcounter+minR,10,4,2,1) );
      jetfindcounter++;
      tempjet[tempjetcounter]->analyze(pythia.event);/*
      if (tempjet[tempjetcounter]->sizeJet()>0&&myJets.size()>0)
      {
        jets.resize(tempjet[tempjetcounter]->sizeJet());
        for(unsigned i=0; i<jets.size(); i++){
          jets[i].pT = tempjet[tempjetcounter]->pT(i);
          jets[i].y= tempjet[tempjetcounter]->y(i);
          jets[i].phi= tempjet[tempjetcounter]->phi(i);
          jets[i].mult= tempjet[tempjetcounter]->multiplicity(i);
          jets[i].r = step*jetfindcounter+minR;
        }
        //printJets(fats,"fats");*/
        myJets = fillJets(myJets,jets);
        //printJets(myJets,"mjets2");
        eventRatio = findRatio(myJets,fats);
     //}
      std::cout<<"Event ratio: "<<eventRatio<<"\n";
    }

    jetFilter(&myJets,fitNUM,fitMAX);
   	tempjets=myJets;
//0
    ipt1=jetMax1(myJets,&fjets[0],0); //1
    fR[0] = myJets[ipt1].r;
    if(fR[0]>0.051)
      std::cout<<"radius: "<<fR[0]<<"\n";
    ipt1=jetMax1(myJets,&fjets[1],fjets[0]);
    fR[1] = myJets[ipt1].r;
//1
    
    for(unsigned i=0; i<myJets.size();++i){
      myJets[i].pT=myJets[i].pT-randomPositive(15,10)/(1-TMath::Power((1-((myJets[i].r*myJets[i].r)/(1.1*1.1))),.5));
     }
    ipt1 = jetMax1(myJets, &fjets[2],0);
    fR[2] = myJets[ipt1].r;
    ipt1= jetMax1(myJets, &fjets[3],fjets[2]);
    fR[3] = myJets[ipt1].r;
       //eventType[iAlag++] = eType(&p1,&p2,myJets[ipt1].phi,myJets[ipt1].y);
    myJets=tempjets;
//2
    
    for(unsigned i=0; i<myJets.size();++i){
       myJets[i].pT=myJets[i].pT-randomPositive(15,10)*myJets[i].r*20*myJets[i].r;
     }
    ipt1 = jetMax1(myJets, &fjets[4],0);
    fR[4] = myJets[ipt1].r;
    ipt1= jetMax1(myJets, &fjets[5],fjets[4]);
    fR[5] = myJets[ipt1].r;
       //eventType[iAlag++] = eType(&p1,&p2,myJets[ipt1].phi,myJets[ipt1].y);
    myJets=tempjets;
//out
      
    for(int i=0; i<preN; i++){
      Xjs[i].Xj = fixXj(fjets[2*i+1],fjets[2*i]); // change Xjs to vector
      //Xjs[i].type = eventType[i];
      Xjs[i].r = maxFloat(fR[2*i+1],fR[2*i]);
      
    }
    /*eventRadius = maxFloat(fR[0],fR[1]);  // does the by radius data colelction nR may be buggy so not using until needed
    for(int i=preN; i<nR;i++){
    	Xjs[i].r =minR+step*(i-preN-1);
    	Xjs[i].Xj=0;
    	if(Xjs[i].r==eventRadius)
    		Xjs[i].Xj=Xjs[0].Xj;
    }*/
   /* switch(Xjs[1].type){
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
    }*/

    t->Fill();
  }
  t->Write();
  f->Write();
  f->Close();
  delete f;
  f=NULL;
}

int main(int argc, char *argv[]){
	std::string filename;
	int fitNUM, fitMAX;
	bool lowpT;
	int nEvent = 1000;
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
