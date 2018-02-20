
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
  float pT =-1;
  float phi =7;
  float y;
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
bool jetEquality(Jet j1, Jet j2);
void printJetR(std::vector<Jet> v, std::string title);

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
float fixXj(Jet j1, Jet j2){
  if(j1.pT>j2.pT)
    return j2.pT/j1.pT;
  else
    return j1.pT/j2.pT;
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
/* finds the highest pT jet from a group of jets with pT less than max*/
void jetMax1(std::vector<Jet> jets, Jet* rJet){ // strong chance this is broken
  assert(jets.size()>0);
  float max=0;
  for(unsigned i=0; i<jets.size(); i++){
    if (jets[i].pT>max)
    {
      max=jets[i].pT;
      rJet->pT= max;
      rJet->r= jets[i].r;
      rJet->phi= jets[i].phi;
      rJet->y = jets[i].y;
      rJet->fatratio = jets[i].fatratio;
      rJet->mult= jets[i].mult;
    }
  }
}
void jetMax1(std::vector<Jet> jets, Jet* rJet, float top){ // strong chance this is broken
  assert(jets.size()>0);
  float max=0;
  for(unsigned i=0; i<jets.size(); i++){
    if (jets[i].pT>max&&jets[i].pT!=top)
    {
      max=jets[i].pT;
      rJet->pT= max;
      rJet->r= jets[i].r;
      rJet->phi= jets[i].phi;
      rJet->y = jets[i].y;
      rJet->fatratio = jets[i].fatratio;
      rJet->mult= jets[i].mult;
    }
  }
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
  assert(fat.size()>1);
  std::vector<Jet> out(jet.size());
  std::vector<float> dR(fat.size());
  int mindex =-1;
  for(unsigned i=0; i<jet.size();i++){  
    for(unsigned j=0; j<fat.size();j++){  // find the differnce in radius between the current jet and all the larger jets
      dR[j] = delR(jet[i],fat[j]);
    } 
    mindex = minimum(dR); // locate the large jet in the same direction as the current jet
    jet[i].fatratio = jet[i].pT/fat[mindex].pT; // set their fat ratio
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
  for(unsigned i=0; i<tempout.size(); i++){
    if(minR<0||tempout[i].fatratio<minR)
      minR = tempout[i].fatratio;
  }
	return minR;
}

std::vector<Jet> fillJets(std::vector<Jet> myJets, std::vector<Jet> Tjets){
  std::vector<Jet> out;
  if(!nullJets(myJets)){
    std::vector<float> dR(Tjets.size());
    int mindex;
    bool contains=false;
    for(unsigned i=0; i<myJets.size(); i++){
      if(myJets[i].fatratio<.9){ //if jet does not have enough energy
      // std::cout<<"delR: ";
        for(unsigned j=0; j<Tjets.size();j++){  //find corresponding higher jet
          dR[j] = delR(myJets[i],Tjets[j]);  //matching by dR but maybe match by within cone is better ?
        // std::cout<<dR[j]<<" ";
         //std::cout<<"loop if for"<<std::endl;
        }
        mindex = minimum(dR);
        unsigned myjetdupcount=0;
        contains=false;
        while(!contains&&myjetdupcount<myJets.size()){
          contains=jetEquality(Tjets[mindex],myJets[myjetdupcount]); //some of these jets need to go to the out
          myjetdupcount++;
        }
        if (contains)
        {
          myjetdupcount--;
          Tjets[mindex].r=myJets[myjetdupcount].r;
          contains=false;
        }
        myjetdupcount=0; 
        while(!contains&&myjetdupcount<out.size()){
          contains=jetEquality(Tjets[mindex],out[myjetdupcount]);
          myjetdupcount++;
        }
        if(!contains){
          out.push_back(Tjets[mindex]);
        }
      //std::cout<<"\n";
      }
      else{
        out.push_back(myJets[i]); // this else might be bad now
        //need to make sure that the other jets we find are not in this jet's Fat jet
      }
    }
  }
  else{
    if(!nullJets(Tjets))
      out=Tjets;
    else{
      cout<<"invalid input in function filljets"<<std::endl;
      throw std::invalid_argument("too many nullJets");
    }
  }
  printJetR(out,"outR");
  printJets(out,"out");
  return out;
}

bool jetEquality(Jet j1, Jet j2){ //really should have used classes 
  assert(j1.pT>=0); //dont check unitialized values
  assert(j2.pT>=0);
  assert(j1.phi<7);
  assert(j2.phi<7);
  bool r=false;
 // cout<<"eq pt1: "<<j1.pT<<" pt2: "<<j2.pT<<std::endl;
  if (j1.pT==j2.pT&&j1.phi==j2.phi){
    //cout<<"jets equal"<<std::endl;
    r=true;
  }
  return r;
}

std::vector<Jet> removeDuplicate(std::vector<Jet> in){
  std::vector<Jet> v(0);
  bool contains=false;
  for(unsigned i=0; i<in.size();i++){
    contains=false;
    for (unsigned j = 0; j < v.size(); ++j){
      contains = jetEquality(in[i],v[j]);
      if(contains){
        break;
      }
    }
    if(!contains){
      v.push_back(in[i]);
    }
  }
  return v;
}

std::vector<Jet> jetFinalFilter(std::vector<Jet> jets, int fitNUM, int fitMAX){ 
  Jet temp;
  jetMax1(jets,&temp,0);
  bool pass= temp.pT>fitNUM&&temp.pT<fitMAX;
  std::vector<Jet> v(0);
  if(pass){
    for (unsigned i = 0; i < jets.size(); ++i)
    {
      if (jets[i].pT>=25)
      {
        v.push_back(jets[i]);
      }
    }
  }
  return v;
}
std::vector<Jet> fatjetFilter(std::vector<Jet> jets, int fitNUM){ 
  bool pass= jets[0].pT>=fitNUM;
  if(pass){
    return jets;
  }
  else{
    std::vector<Jet> v(0);
    return v;
  }
}
std::vector<Jet> jetFilter(std::vector<Jet> jets, int fitNUM, int fitMAX){ 
  Jet temp;
  jetMax1(jets,&temp,0);
  bool pass= temp.pT>fitNUM&&temp.pT<fitMAX;
  if(!pass){
    std::vector<Jet> no(0);
    return no;
  }
  else{
    return jets;
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
inline void printJetR(std::vector<Jet> v, std::string title){
  std::cout<<title<<" Rs: ";
  for(unsigned i=0; i<v.size();i++){
    std::cout<<v[i].r<<" ";
  }
  std::cout<<std::endl;
}
void clear(std::vector<Jet>* v){
  while(!v->empty()){
    v->pop_back();
  }
}

inline float maxJet(std::vector<Jet> v){
  float r=0;
  for (unsigned i = 0; i < v.size(); ++i)
  {
    if(r<v[i].pT)
      r=v[i].pT;
  }
  return r;
}

void fillXjs(XjT* xT, std::vector<Jet> myJets){
  printJets(myJets,"filling Xj");
  assert(myJets.size()>1);
  Jet jettemp1;
  jetMax1(myJets,&jettemp1,0);
  assert(jettemp1.pT>0);
  cout<<&jettemp1<<": "<<jettemp1.pT;
  Jet jettemp2;
  jetMax1(myJets,&jettemp2,jettemp1.pT);
  assert(jettemp2.pT>0);
  cout<<" "<<&jettemp2<<": "<<jettemp2.pT<<" Xj: ";
  xT->Xj = fixXj(jettemp1,jettemp2);
  cout<<xT->Xj<<"\n";
  xT->r = maxFloat(jettemp1.r,jettemp2.r);
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
  const float maxR = .4;
  const float maxloop = (maxR-minR)/step;
  std::vector<SlowJet*> tempjet(0);
  SlowJet* deltemp = NULL;
  int tempjetcounter;
  Jet jettemp1;
  const int nXj = 3;
  XjT Xjs[nXj];
  t->Branch("Xj", &Xjs[0].Xj); 
  t->Branch("XjQ0", &Xjs[1].Xj);
  t->Branch("XjQ1", &Xjs[2].Xj);
  t->Branch("LeadR0",&Xjs[0].r);
  float eventRatio=0;
  short jetfindcounter=0;
  std::vector<Jet> jets;
  std::vector<Jet> myJets;
  std::vector<Jet> fats;
  std::vector<Jet> tempjets;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
  	if(iEvent%30==0)  
  		cout<<"Event N: "<<iEvent<<'\n';
    if (!pythia.next()){
      cout<<"pythia.next() failed"<<"\n";
      nEvent--;
      continue;
    }
    //std::cout<<"event loop"<<iEvent<<'\n';

    clear(&myJets);
    clear(&fats);
    clear(&jets);
    tempjet.push_back( new SlowJet(-1,1,10,4,2,1));  //set up the comparison jets
    tempjet[0]->analyze(pythia.event);
    for(int i=0; i<tempjet[0]->sizeJet(); i++){
      jettemp1.pT = (float) tempjet[0]->pT(i);
      jettemp1.y= (float) tempjet[0]->y(i);
      jettemp1.phi= (float) tempjet[0]->phi(i);
      jettemp1.mult= (float) tempjet[0]->multiplicity(i);
      jettemp1.r = 1;
      fats.push_back(jettemp1);
    }
    printJets(fats,"unfiltered fat");
    fats=fatjetFilter(fats,fitNUM);
    printJets(fats, "filtered fat");
    if(fats.size()<=1){
      nEvent--;
      continue;
    }

    tempjet.push_back(new SlowJet(-1,minR,10,4,2,1)); // set up the smallest jets 
    jetfindcounter=1;
    tempjet[1]->analyze(pythia.event);
    for(int i=0; i<tempjet[1]->sizeJet(); i++){
      jettemp1.pT = (float)tempjet[1]->pT(i);
      jettemp1.y= (float)tempjet[1]->y(i);
      jettemp1.phi= (float)tempjet[1]->phi(i);
      jettemp1.mult= (float)tempjet[1]->multiplicity(i);
      jettemp1.r = minR;
      myJets.push_back(jettemp1);
    }

    eventRatio=0;
    tempjetcounter=1;
    while(eventRatio<.9&&jetfindcounter<maxloop){ 
      tempjetcounter++;
      tempjet.push_back(new SlowJet(-1,step*jetfindcounter+minR,10,4,2,1) );
      jetfindcounter++;
      assert(tempjetcounter<tempjet.size());
      assert(tempjet[tempjetcounter]!=NULL);
      tempjet[tempjetcounter]->analyze(pythia.event);
      /*for(int i=0; i<tempjet.size(); i++){
        cout<<tempjet[i]<<": "<<tempjet[i]->sizeJet()<<"\n";
      }*/
      if (tempjet[tempjetcounter]->sizeJet()>1)
      {
        clear(&jets);
        for(int i=0; i<tempjet[tempjetcounter]->sizeJet(); i++){
          jettemp1.pT = (float) tempjet[tempjetcounter]->pT(i);
          jettemp1.y= (float) tempjet[tempjetcounter]->y(i);
          jettemp1.phi= (float) tempjet[tempjetcounter]->phi(i);
          jettemp1.mult= (float) tempjet[tempjetcounter]->multiplicity(i);
          jettemp1.r = (float) step*jetfindcounter+minR;
          jets.push_back(jettemp1);
        }
        if(myJets.size()>1){
          printJets(myJets,"prefill");
          printJets(jets,"new");
          myJets = fillJets(myJets,jets);
          printJets(myJets, "post fill");
          printJets(fats,"fats");
          /*
          cout<<"myJets: \n";
          for(unsigned i=0; i<myJets.size(); i++){
          cout<<&myJets[i]<<": "<<myJets[i].pT<<", "<<myJets[i].r<<"\n";
          }
          cout<<"new jets: \n";
          for(unsigned i=0; i<jets.size(); i++){
            cout<<&jets[i]<<": "<<jets[i].pT<<", "<<jets[i].r<<"\n";
          }*/
        //printJets(myJets,"mjets2");
        }
        else{
          myJets=jets;
        }
        eventRatio = findRatio(myJets,fats);
      }
      else{
        deltemp = tempjet.back();
        tempjet.pop_back();
        delete deltemp;
        deltemp=NULL;
        tempjetcounter--;
      }
      std::cout<<"Event ratio: "<<eventRatio<<"\n";
    }
    while(!tempjet.empty()){
      deltemp = tempjet.back();
      tempjet.pop_back();
      delete deltemp;
      deltemp=NULL;
    }
    printJets(myJets,"before filter");
    myJets =jetFinalFilter(myJets,fitNUM,fitMAX);
    if(myJets.size()<=1){
      nEvent--;
      continue;
    }
    printJets(myJets,"after filter");
   	tempjets=myJets;
//0
    fillXjs(&Xjs[0],myJets);
//1
    float tempenergy;
    
    for(unsigned i=0; i<myJets.size();++i){ 
      myJets[i].pT=myJets[i].pT * TMath::Power(E,-1*randomPositive(2,1)*myJets[i].r*myJets[i].r);
     }
    printJets(myJets,"fill 2");
    fillXjs(&Xjs[1],myJets);
    myJets=tempjets;
//2   cuts are making negative numbers
    
    for(unsigned i=0; i<myJets.size();++i){
       myJets[i].pT=myJets[i].pT-randomPositive(15,10)*myJets[i].r*myJets[i].r;
       cout<<"3r: "<<myJets[i].r<<" ";
    }
    printJets(myJets,"fill 3");
    fillXjs(&Xjs[2],myJets);
    myJets=tempjets;

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