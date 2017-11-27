// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details
/*
Eventually
 time to start looking at b quark and path length
 */
/* There are currently 6 jet quenching algorithms that select for parton called 1,2,3,A,B,C
1 = -GAUS(0,20)
2= -GAUS(10,20)
3= -GAUS(20,20)
A = *(1-GAUS(.2,.05))
B= *(1-GAUS(.2, .1))
C= *(1-GAUS(.2,.2))
 These algorithms do not select for parton (yet) and are
4= -GAUS(20, 30)
5= - GAUS(20,40)
ZR = -GAUS(10,10) - rateB
ZR2 = -GAUS(20,20) - rateB
ZR3 = xGAUS(.2,.3)-rateB
Z4 = xGAUS(.2,.1)-GAUS(20,30);
D = *(1-GAUS(.3,.2))
E= *(1-GAUS(.1,.2))
P = Interpolation of decreasing fractional energy loss wrt pT
   Rate algorithms subtract a set amount of energy per charged particle
rateA = 1.5
rateB = 2;
rateC = 1;
These are agorlithms are displayed as totals and w.r.t. flavor
*/


#include "Pythia8/Pythia.h"
using namespace Pythia8;
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "spline.h"
#include <fstream>

int mycount=0;
Double_t E= 2.71828182845904523536;
struct MJet{
  float pT, phi, y;
  int evntN, jetN;
  vector<int> partN;
};

struct Parton{
  float px, py, pz, eta;
  int id;
};
struct XjT{
  float Xj;
  float fat;
  int type;
};

float randomPositive(double mean, double sig){
  TRandom3 *tr3= new TRandom3(0);
  float r = tr3->Gaus(mean, sig);
  if(r>0){
    return r;
  }
  else{
    return randomPositive(mean, sig);
  }
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
float delPhi(float phi1, Parton p1){
  return TMath::Abs(phi1 - TMath::ATan2(p1.py,p1.px));
}

float quadrature(float x, float y , float z){
  return TMath::Sqrt(x*x+y*y+z*z);
}


float etaToRad(float eta1){
  Double_t e = (double) eta1*-1;
  return  2* TMath::ATan((double)TMath::Power(E,e));
}

float delR(float phi1, float y1, Parton p1){
  y1 = etaToRad(y1);
  float y2 = etaToRad(p1.eta);
  phi1 = delPhi(phi1, p1);
  return quadrature(TMath::Abs(y1-y2),phi1,0);
}

void swapPart(Parton *p1, Parton *p2){
  Parton temp = *p1;
  *p1 = *p2;
  *p2=temp;
}

void setLead(Parton *p1, Parton *p2, float phi, float y){
  if(delR(phi,y,*p1)>delR(phi,y,*p2))
    {
      swapPart(p1, p2);

    }
}
bool isQuark(Parton p1){
  if(TMath::Abs(p1.id)>=1&& TMath::Abs(p1.id)<=8)
    return true;
  return false;
}

int eType(Parton *p1,Parton *p2, float phi, float y){
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


int jetMax1(std::vector<float> jets, float* jet, float max){
  *jet=0;
  int r=-1;
  if(max==0){
      for(int i=0; i<jets.size();i++){
	       if(jets[i]>*jet){
	       *jet=jets[i];
         r=i;
	       }
      }
  }
  else{
    for(int i=0; i<jets.size();i++){
	     if(jets[i] >*jet && jets[i]<max){
	       *jet=jets[i];
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

void addJet(float *jet1, float *phi1, float jet2, float phi2){
  float x = (*jet1)* TMath::Cos((*phi1)) + jet2 * TMath::Cos(phi2);
  float y = (*jet1)*TMath::Sin((*phi1)) + jet2 * TMath::Sin(phi2);
  *jet1 = TMath::Sqrt((double)(x*x+y*y));
  *phi1 = TMath::ATan2(y,x);
}

void addJet2(float *jet1, float phi1, float *jet2, float phi2, float jet3, float phi3){
  if(TMath::Abs(phi1-phi3)<TMath::Abs(phi2-phi3)){ //if delphi for J1-J3 is less than delphi J2-J3
    addJet(jet1, &phi1,jet3, phi3);

  }
  else{
    addJet(jet2, &phi2, jet3,phi3);
  }
}

void printP(vector<int> parts, vector<int> particles){
  for(unsigned int i=0; i<parts.size(); i++){
    cout<<particles[parts[i]]<<endl;
  }
}

void printParton(Parton p){
  cout<<p.id<<": Px: "<<p.px<<" Py: "<<p.py<<" Pz: "<<p.pz;
}

void setJetData(int index, float* phi, float* y, float jet_phi[], float jet_y[]){
  if(index==-1){
    *phi=0;
    *y=0;
  }
  *phi = jet_phi[index];
  *y = jet_y[index];
}
/**
arg#0 = run command
arg#1 = run state low for lowpT high for highpT
arg#2 = run state o for over write a for append 
arg#3 = number of events 
can be run without arguments 
**/
int main(int argc, char *argv[]){ 
  int nEvent =2000;
  bool lowpT =true;
  int fitNUM, fitMAX;
  std::string outfilebegin;
  std::string outfileend = ".root";
  std::string writeSet = "RECREATE";
  int fileN =1;
  if(argc==2||argc==3){
    std::string arg1(argv[1]);
    if(arg1=="low")
      lowpT=true;
    else if(arg1=="high")
      lowpT=false;
    else
      return 8;
  }
  if(lowpT)
    outfilebegin = "dijet100";
  else
    outfilebegin = "dijet200";

  if(argc<3||(argc==3&&argv[2][0]=='o')){
    bool deleting=true;
    using namespace std;
    string filename;
    while(deleting){
      filename=outfilebegin+to_string(fileN)+outfileend;
      if(remove(filename.c_str())!=0){
        deleting=false;
      }
      else{
       fileN++; 
      }
    }
    fileN=1;
  }
  else{
    bool running=true;
    using namespace std;
    string filename;
    while(running){
      filename=outfilebegin+to_string(fileN)+outfileend;
      ifstream infile(filename);
      if(infile.good()){
        fileN++;
      }
      else{
        running=false;
      }
    }
  }
  std::string outfilename = outfilebegin+std::to_string(fileN)+outfileend;
  TFile *f = new TFile(outfilename.c_str(),writeSet.c_str());
  TTree *t=new TTree("tree100","100pThat events");
  if(lowpT){
    fitNUM =100;
    fitMAX = 126;
  }
  else{
    fitNUM=200;
    fitMAX=3000;
  }
  //Interpolation values
  
  std::vector<double> interE = {50,70,90,100,110,126,140,170,200}; 
  std::vector<double> interM={.5,.43,.35,.3,.25,.2,.15,.1,.05};
  tk::spline spline;
  spline.set_points(interE,interM);

  Pythia pythia;
  pythia.readString("Beams:eCM = 2760.");
  pythia.readString("HardQCD:all = on");
  if(lowpT)
    pythia.readString("PhaseSpace:pTHatMin = 75.");
  else
    pythia.readString("PhaseSpace:pTHatMin = 150.");
  pythia.init();

  SlowJet fatjet(-1, 0.8, 10,4,2,1);
  SlowJet slowJet( -1, 0.4, 10, 4, 2, 1);

  int jet_n;
  int fatjetcount;
  float highjet;
  float lowjet;
  float jet_pt[100]; // do I need fatjets?

  float jet_y[100];
  float jet_phi[100];
  float mult[100];
  float fatmult[100];
  Parton p1; Parton p2;
  float  Xj,X1,X2,X3,XA,XB,QQ1,QQ2,QQ3,QQA,QQB,QQC,QG1,QG2,QG3,QGA,QGB,QGC,GQ1,GQ2,GQ3,GQA,GQB,GQC,GG1,GG2,GG3,GGA,GGB, GGC, XC, ZR, ZR2,ZR3,Z4,RAQQ, RAQG, RAGQ, RAGG, RBQQ,RBQG,RBGQ,RBGG,RCQQ,RCQG,RCGQ,RCGG;
  float xrate, xrateB, xrateC;
  float X4,X5,XD,XE,XP;

  const int nXj = 20;
  XjT Xjs[nXj];

  t->Branch("highjet",&highjet);t->Branch("lowjet",&lowjet);
  t->Branch("Xj", &Xj);t->Branch("X1", &X1);t->Branch("X2", &X2); t->Branch("X3", &X3);t->Branch("XA", &XA); t->Branch("XB", &XB);t->Branch("XC", &XC);
  t->Branch("X4", &X4);t->Branch("X5", &X5);
  t->Branch("QQ1",&QQ1);t->Branch("QQ2",&QQ2);t->Branch("QQ3",&QQ3);t->Branch("QQA",&QQA);t->Branch("QQB", &QQB);t->Branch("QQC",&QQC);t->Branch("QG1",&QG1);t->Branch("QG2",&QG2);t->Branch("QG3",&QG3);t->Branch("QGA",&QGA);t->Branch("QGB",&QGB);t->Branch("QGC",&QGC);t->Branch("GQ1",&GQ1);t->Branch("GQ2",&GQ2);t->Branch("GQ3",&GQ3);t->Branch("GQA",&GQA);t->Branch("GQB",&GQB);t->Branch("GQC",&GQC);t->Branch("GG1",&GG1);t->Branch("GG2",&GG2);t->Branch("GG3",&GG3);t->Branch("GGA",&GGA);t->Branch("GGB",&GGB);t->Branch("GGC", &GGC);
  t->Branch("XD",&XD);t->Branch("XE", &XE);t->Branch("XP",&XP);
  t->Branch("xrate", &xrate);t->Branch("xrateB",&xrateB);t->Branch("xrateC", &xrateC);
  t->Branch("ZR", &ZR);t->Branch("ZR2", &ZR2);t->Branch("ZR3",&ZR3);t->Branch("Z4",&Z4);
  t->Branch("RAQQ", &RAQQ);t->Branch("RBQQ",&RBQQ); t->Branch("RCQQ",&RCQQ);t->Branch("RAQG",&RAQG);t->Branch("RBQG",&RBQG);t->Branch("RBGQ",&RBGQ);t->Branch("RBGG",&RBGG);t->Branch("RCQG", &RCQG); t->Branch("RAGG",&RAGG);t->Branch("RAGQ",&RAGQ);t->Branch("RAGG",&RAGG);t->Branch("RBQG",&RBQG);t->Branch("RCGQ",&RCGQ);t->Branch("RBGG", &RBGG); t->Branch("RCGG",&RCGG);
  std::vector<float> jetpT;
  std::vector<float> tempjets;
  std::vector<float> fatTemp;
  std::vector<float> fatpT;
  std::string temp;
  for(int i=0; i<nXj; i++){
    temp = "fat"+std::to_string(i);
    t->Branch(temp.c_str(), &Xjs[i].fat);
  }

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next())
      continue;
    slowJet. analyze( pythia.event );
    fatjet.analyze(pythia.event);
    p1.px = pythia.event[5].px();
    p1.py=pythia.event[5].py();
    p1.eta = pythia.event[5].eta();
    p2.eta = pythia.event[6].eta();
    p2.px=pythia.event[6].px();
    p2.py=pythia.event[6].py();
    p1.id= pythia.event[5].id();
    p2.id= pythia.event[6].id();
    jet_n=0;
    fatjetcount=0;
    jetpT.resize(slowJet.sizeJet());
    tempjets.resize(slowJet.sizeJet());
    fatTemp.resize(fatjet.sizeJet());
    fatpT.resize(fatjet.sizeJet());
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      jet_pt[ jet_n ] = slowJet.pT(i);
      jet_y[ jet_n ] = slowJet.y(i);
      jet_phi[ jet_n ] = slowJet.phi(i);
      mult[jet_n] = slowJet.multiplicity(i);
      tempjets[i] = slowJet.pT(i);
      jet_n++;
    }
    for(int i=0; i<fatjet.sizeJet();++i){ //will I need the phi eta to make parton selections??
      fatTemp[i]=fatpT[fatjetcount]=fatjet.pT(i);
      fatmult[fatjetcount++]=fatjet.multiplicity(i);
    }
    
    jetpT = tempjets;
    const int nfjets=40;
    const float RATE = 1.5;
    const float RATEB=2;
    const float RATEC=1;
    float fjets[nfjets];
    float leadfatjets[nfjets/2];
    int eventType[nfjets/2];
    int iAlag=0;
    Xj=0; QQ1=0;QQ2=0;QQ3=0;QQA=0;QQB=0;QG1=0;QG2=0;QG3=0;QGA=0;QGB=0;GQ1=0;GQ2=0;GQ3=0;GQA=0;GQB=0;GG1=0;GG2=0;GG3=0;GGA=0;GGB=0;xrate=0;XC=0;QQC=0;QGC=0;GQC=0;GGC=0;RAQQ=0; RAQG=0; RAGQ=0; RAGG=0; RBQQ=0;RBQG=0;RBGQ=0;RBGG=0;RCQQ=0;RCQG=0;RCGQ=0;RCGG=0;

     int ipt1=0;
     int ipt2=1;
    
    highjet=fjets[0] = jet_pt[0];
    lowjet =fjets[1] = jet_pt[1];
    leadfatjets[0]=fatpT[0];

    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jet_pt[i]-randomPositive(0,20);
       fatpT[i]=fatpT[i]-randomPositive(0,20);
     }

       ipt1 = jetMax1(jetpT, &fjets[2],0);
       jetMax1(jetpT, &fjets[3],fjets[2]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;

       ipt1 = jetMax1(fatpT, &leadfatjets[1],0);
       fatpT=fatTemp;

    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jet_pt[i]-randomPositive(10,20);
       fatpT[i]=fatpT[i]-randomPositive(10,20);
     }
       ipt1=jetMax1(jetpT, &fjets[4],0);
       jetMax1(jetpT, &fjets[5],fjets[4]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[2],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jet_pt[i]-randomPositive(20,20);
       fatpT[i]=fatpT[i]-randomPositive(20,20);
     }
       ipt1=jetMax1(jetpT, &fjets[6],0);
       jetMax1(jetpT, &fjets[7],fjets[6]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[3],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.2, .05));
       fatpT[i]=fatpT[i]*(1-randomPositive(.2, .05));
     }
       ipt1=jetMax1(jetpT, &fjets[8],0);
       jetMax1(jetpT, &fjets[9],fjets[8]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[4],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.2, .1));
       fatpT[i]=fatpT[i]*(1-randomPositive(.2, .1));
     }
       ipt1=jetMax1(jetpT, &fjets[10],0);
       jetMax1(jetpT, &fjets[11],fjets[10]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[5],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]-(mult[i]*RATE);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATE);
     }
       ipt1=jetMax1(jetpT, &fjets[12],0);
       jetMax1(jetpT, &fjets[13],fjets[12]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[6],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.2, .2));
       fatpT[i]=fatpT[i]*(1-randomPositive(.2, .2));
     }
       ipt1=jetMax1(jetpT, &fjets[14],0);
       jetMax1(jetpT, &fjets[15],fjets[14]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[7],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]-randomPositive(20,30);
       fatpT[i]=fatpT[i]-randomPositive(20,30);
     }
       ipt1=jetMax1(jetpT, &fjets[16],0);
       jetMax1(jetpT, &fjets[17],fjets[16]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[8],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]-randomPositive(20,40);
       fatpT[i]=fatpT[i]-randomPositive(20,40);
     }
       ipt1=jetMax1(jetpT, &fjets[18],0);
       jetMax1(jetpT, &fjets[19],fjets[18]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[9],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.2, .2)); //this is a duplicate 
       fatpT[i]=fatpT[i]*(1-randomPositive(.2, .2));
     }
       ipt1=jetMax1(jetpT, &fjets[20],0);
       jetMax1(jetpT, &fjets[21],fjets[20]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[10],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]-(mult[i]*RATEB);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEB);
     }
       ipt1=jetMax1(jetpT, &fjets[22],0);
       jetMax1(jetpT, &fjets[23],fjets[22]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[11],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]-(mult[i]*RATEC);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEC);
     }
       ipt1=jetMax1(jetpT, &fjets[24],0);
       jetMax1(jetpT, &fjets[25],fjets[24]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[12],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){ // these are for the Z algorithms 
       jetpT[i]=jetpT[i]-(mult[i]*RATEB);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEB);
     }
       ipt1=jetMax1(jetpT, &fjets[26],0);
       jetMax1(jetpT, &fjets[27],fjets[26]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[13],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){ // ditto
       jetpT[i]=jetpT[i]-(mult[i]*RATEB);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEB);
     }
       ipt1=jetMax1(jetpT, &fjets[28],0);
       jetMax1(jetpT, &fjets[29],fjets[28]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[14],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){ //ditto
       jetpT[i]=jetpT[i]-(mult[i]*RATEB);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEB);
     }
       ipt1=jetMax1(jetpT, &fjets[30],0);
       jetMax1(jetpT, &fjets[31],fjets[30]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[15],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){ // ditto
       jetpT[i]=jetpT[i]-(mult[i]*RATEB);
       fatpT[i]=fatpT[i]-(fatmult[i]*RATEB);
     }
       ipt1=jetMax1(jetpT, &fjets[32],0);
       jetMax1(jetpT, &fjets[33],fjets[32]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[16],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.3,.2));
       fatpT[i]=fatpT[i]*(1-randomPositive(.3,.2));
     }
       ipt1=jetMax1(jetpT, &fjets[34],0);
       jetMax1(jetpT, &fjets[35],fjets[34]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[17],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(.1,.2));
       fatpT[i]=fatpT[i]*(1-randomPositive(.1,.2));
     }
       ipt1=jetMax1(jetpT, &fjets[36],0);
       jetMax1(jetpT, &fjets[37],fjets[36]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[18],0);
       fatpT=fatTemp;
    for(int i=0; i<slowJet.sizeJet();++i){
       jetpT[i]=jetpT[i]*(1-randomPositive(spline(jetpT[i]),.2));
       fatpT[i]=fatpT[i]*(1-randomPositive(spline(fatpT[i]),.2));
     }
       ipt1=jetMax1(jetpT, &fjets[38],0);
       jetMax1(jetpT, &fjets[39],fjets[38]);
       eventType[iAlag++] = eType(&p1,&p2,jet_phi[ipt1],jet_y[ipt1]);
       jetpT=tempjets;
       ipt1 = jetMax1(fatpT, &leadfatjets[19],0);
       fatpT=fatTemp;
    std::string fatBranchName;
    for(int i=0; i<nXj; i++){
      Xjs[i].Xj = fixXj(&fjets[2*i+1],&fjets[2*i],fitNUM,fitMAX);
      Xjs[i].type = eventType[i];
      Xjs[i].fat = fjets[2*i] / leadfatjets[i];
      if(Xjs[i].Xj>1){
	       Xjs[i].Xj = 1/Xjs[i].Xj;
	       if(Xjs[i].type ==2)
	           Xjs[i].type=3;
	       else if(Xjs[i].type==3)
	           Xjs[i].type=2;
      }
    }
    Xj = Xjs[0].Xj;
    X1 = Xjs[1].Xj;
    X2 = Xjs[2].Xj;
    X3 = Xjs[3].Xj;
    XA = Xjs[4].Xj;
    XB = Xjs[5].Xj;
    xrate = Xjs[6].Xj;
    XC = Xjs[7].Xj;
    X4 = Xjs[8].Xj;
    X5 = Xjs[9].Xj;
    xrateB = Xjs[11].Xj;
    xrateC = Xjs[12].Xj;
    ZR = Xjs[13].Xj;
    ZR2 = Xjs[14].Xj;
    ZR3 = Xjs[15].Xj;
    Z4 = Xjs[16].Xj;
    XD = Xjs[17].Xj;
    XE = Xjs[18].Xj;
    XP = Xjs[19].Xj;
    switch(Xjs[1].type){
    case 1:
	     QQ1 = Xjs[1].Xj;
       QQ2=Xjs[2].Xj;
	     QQ3=Xjs[3].Xj;
	     QQA=Xjs[4].Xj;
	     QQB=Xjs[5].Xj;
	     QQC=Xjs[10].Xj;
	     RAQQ = xrate;
	     RBQQ = Xjs[11].Xj;
	     RCQQ = Xjs[12].Xj;
	     break;
    case 2:
       QG1 = Xjs[1].Xj;
	     QG2=Xjs[2].Xj;
	     QG3=Xjs[3].Xj;
	     QGA=Xjs[4].Xj;
	     QGB=Xjs[5].Xj;
	     QGC=Xjs[10].Xj;
	     RAQG = xrate;
	     RBQG = Xjs[11].Xj;
	     RCQG = Xjs[12].Xj;
	     break;
    case 3:
       GQ1 = Xjs[1].Xj;
	     GQ2=Xjs[2].Xj;
	     GQ3=Xjs[3].Xj;
	     GQA=Xjs[4].Xj;
	     GQB=Xjs[5].Xj;
	     GQC=Xjs[10].Xj;
	     RAGQ = xrate;
	     RBGQ = Xjs[11].Xj;
	     RCGQ = Xjs[12].Xj;
	     break;
    case 4:
       GG1 = Xjs[1].Xj;
	     GG2=Xjs[2].Xj;
	     GG3=Xjs[3].Xj;
	     GGA=Xjs[4].Xj;
	     GGB=Xjs[5].Xj;
	     GGC=Xjs[10].Xj;
	     RAGG = xrate;
	     RBGG = Xjs[11].Xj;
	     RCGG = Xjs[12].Xj;
	     break;
    }
    t->Fill();
  }
  t->Write();
  f->Write();
  f->Close();
  //cout<<mycount<<std::endl; //for debug
  return 0;
}