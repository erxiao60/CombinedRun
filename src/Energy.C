#define Energy_cxx
#include "RootConvert.h"
#include <TStyle.h>
#include <sstream>
#include <string.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <stdio.h>
#include <TCanvas.h>
#include <TSpectrum.h>

#include "TFile.h"
#include <TF1.h>
#include "TInterpreter.h"
#include <TStyle.h>
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TProfile.h"
using namespace std;
using namespace cosmictest;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TrackSeek(Int_t nHits,Int_t *Layer,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][2],Int_t *Dy,Int_t *Side,Int_t *GID,Int_t FitTrack[15])
{
  Float_t MaxADC[14];
  Int_t   MaxBar[14];//[2]both sides of BGO bar
  memset(MaxADC,0,sizeof(MaxADC));
  memset(MaxBar,0,sizeof(MaxBar));
  Float_t ADCCutPed=0; 
  for(int i=0; i<nHits;i++)//only check side 0
  {
    if(Dy[i]==8&&Side[i]==0)//only use Dy8 hits to get a track!!!
    {
    if (ADC[i]>(PedestalPar[GID[i]][0]+3*PedestalPar[GID[i]][1]))
    {
      ADCCutPed=ADC[i]-PedestalPar[GID[i]][0];
      if(ADCCutPed>MaxADC[Layer[i]])  
      {
        MaxADC[Layer[i]]=ADCCutPed;
        MaxBar[Layer[i]]=Bar[i];//0--21
      }
    }
    }
  } //single event Loop end
FitTrack[14]=0;

if((MaxADC[0]>200||MaxADC[1]>200)&&(MaxADC[12]>200||MaxADC[13]>200))
{
FitTrack[14]=1;
for(int il=0;il<Nlayer;il++)
{
FitTrack[il]=MaxBar[il];
}}
}

void RootConvert::Loop(char *RootFileName)
{
gROOT->ProcessLine(".L /home/zhzhy/DRAW/bes3plotstyle.C");
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//output .root .eps files
string namesuffix="_Energy.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[50]="./Energy/";
char *EPSFile=strcat(FilePath,outputfile.c_str());
cout<<"Output file: "<<EPSFile<<endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read pedestal data
cout<<"Read Pedestal..."<<endl;

const int nPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
//const int nGID=(Nlayer*Nside)*Nbar*Ndy;
Double_t PedestalPar[nGID][nPar];
ifstream PedPar;
PedPar.open("./Pedestal/Pedestal.txt");
if (!PedPar.good())
{cout<<"Can not open Pedestal TXTFlie  File!!!"<<endl;
exit(-1);
}
for(int ig=0;ig<nGID;ig++)
{
PedPar>>PedestalPar[ig][0];
PedPar>>PedestalPar[ig][1];
}
PedPar.close();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read MIPs data
const int nGbar=Nlayer*Nbar;
Double_t MIPsPar[nGbar*Nside][4];
ifstream M_Par;
M_Par.open("./MIPsV2/MIPsV2.txt");
if(!M_Par.good())
{
cout<<"can not open out put fitted par file :MIPs_LangauPar.txt"<<endl;
exit(1);
}
for(int i=0;i<nGbar*Nside;i++)
for(int j=0;j<4;j++)
M_Par>>MIPsPar[i][j];

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define 
//TH*F
TH1F *L_Energy[Nlayer][3];//[0],[1]:side0,1; 2:square root of side0*side1
TH1F *T_Energy[3];

for (int il=0;il<Nlayer;il++){
  for(int is=0; is<Nside+1;is++){
      // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
    char cc[30];
    if(is==2){
    sprintf(cc,"Layer%d_corrected_energy",il);
    }
    else
    sprintf(cc,"Layer%d_Side%d_Energy",il,is);
    L_Energy[il][is]=new TH1F(cc,cc,50,0,10);
    }
}
  for(int is=0; is<Nside+1;is++){
    char cc[30];
    if(is==2){
    sprintf(cc,"corrected_energy");
    }
    else
    sprintf(cc,"Side%d_Energy",is);
    T_Energy[is]=new TH1F(cc,cc,175,0,350);
    }


Double_t energy[Nlayer][Nbar][3];
for (int il=0;il<Nlayer;il++)
  for (int ib=0;ib<Nbar;ib++)
    for(int is=0; is<Nside+1;is++)
    energy[il][ib][is]=0.;

TF1 *MIPs=new TF1("MIPs","landau",0,350);
TF1 *myMIPs=new TF1("myMIPs",langaufun,0,350,4);
myMIPs->SetParNames("Width","MP","Area","GSigma");
gStyle->SetOptStat(111);
gStyle->SetOptFit(1111);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//LOOP START
cout<<"Loop root file ..."<<endl;
Int_t BadTrack=0;
Int_t Tevent=0;
if (fChain == 0) {cout<<"fChain=0,return!";return;}
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) {cout<<"ientry="<<ientry<<",break";break;}
  nb = fChain->GetEntry(jentry);   
  if(jentry%1000==0)
  cout<<"Total Events:"<<jentry<<endl;
//if(jentry==10000) break;
 Int_t FitTrack[15];//
 Double_t Par[2];
 TrackSeek(nHits,Layer,Bar,ADC,PedestalPar,Dy,Side,GID,FitTrack);
if(FitTrack[14]==1){

  //    cout<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //  cout<<"Tevent:"<<Tevent<<endl;
  for(int i=0;i<nHits;i++){
    if(Dy[i]==8&&Bar[i]!=0&&Bar[i]!=23){//both sides
      if (ADC[i]>(PedestalPar[GID[i]][0]+2*PedestalPar[GID[i]][1])&&ADC[i]<=13000){
       int iGbar=(Layer[i]*Nside+Side[i])*Nbar+Bar[i];
       energy[Layer[i]][Bar[i]][Side[i]]=(ADC[i]-PedestalPar[GID[i]][0])/MIPsPar[iGbar][1];
   //   cout<<"Tevent:"<<Tevent<<"\n"<<"Layer:"<<Layer[i]<<" Bar:"<<Bar[i]<<" Side:"<<Side[i]<<" ADC:"<<ADC[i]<<" MIPs:"<<MIPsPar[iGbar][1]<<" Pedestal:"<<PedestalPar[GID[i]][0]<<endl;
       
      }
      else if(ADC[i]>13000)
      cout<<"Tevent:"<<Tevent<<"\n"<<"Layer:"<<Layer[i]<<" Bar:"<<Bar[i]<<" Side:"<<Side[i]<<" ADC:"<<ADC[i]<<endl;
     }
     }
   //for energy correction
   for (int il=0;il<Nlayer;il++){
     for (int ib=0;ib<Nbar;ib++){
      energy[il][ib][2]=TMath::Sqrt(energy[il][ib][0]*energy[il][ib][1]);
      
      }
    }
   //energy fill
    Double_t Et[3]={0.,0.,0.};
    Double_t El[Nlayer][3];
    for(int is=0;is<Nside+1;is++){
      Et[is]=0.;
      for(int il=0;il<Nlayer;il++){
        El[il][is]=0.;  
        }
     }
    for(int is=0;is<Nside+1;is++){
      for(int il=0;il<Nlayer;il++){
        for(int ib=0;ib<Nbar;ib++){
        if(energy[il][ib][is]>0.01)
        El[il][is]+=energy[il][ib][is];
        }
        Et[is]+=El[il][is];
        L_Energy[il][is]->Fill(El[il][is]);

      }
      T_Energy[is]->Fill(Et[is]);
    }

for (int il=0;il<Nlayer;il++)
  for (int ib=0;ib<Nbar;ib++)
    for(int is=0; is<Nside+1;is++)
    energy[il][ib][is]=0.;
  }

else 
  {BadTrack++;
  }
Tevent++;
}
//Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Double_t MPar[3];
// Energy Fit and Draw
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

TCanvas *EMCEnergy[3];
for(int is=0;is<Nside+1;is++)
{
char d[50];
sprintf(d,"Side%d_Energy",is);
EMCEnergy[is]=new TCanvas(d,d,0,0,600,800);
EMCEnergy[is]->Divide(4,4);
for(int il=0;il<Nlayer;il++)
{
EMCEnergy[is]->cd(il+1);
L_Energy[il][is]->Draw();
MIPs->SetRange(0,5);
L_Energy[il][is]->Fit(MIPs,"Q0");
MIPs->GetParameters(MPar);
myMIPs->SetRange(MPar[1]-3.6*MPar[2],MPar[1]+6.50*MPar[2]);
myMIPs->SetParameter(0,MPar[2]);
myMIPs->SetParameter(1,MPar[1]);
myMIPs->SetParameter(3,MPar[2]);
myMIPs->SetLineColor(2);
myMIPs->SetLineWidth(2);
L_Energy[il][is]->Fit(myMIPs,"R");

}
EMCEnergy[is]->cd(15);
T_Energy[is]->Draw();
MIPs->SetRange(0,35);
T_Energy[is]->Fit(MIPs,"Q0");
MIPs->GetParameters(MPar);
myMIPs->SetRange(MPar[1]-3.6*MPar[2],MPar[1]+6.50*MPar[2]);
myMIPs->SetParameter(0,MPar[2]);
myMIPs->SetParameter(1,MPar[1]);
myMIPs->SetParameter(3,MPar[2]);
myMIPs->SetLineColor(2);
myMIPs->SetLineWidth(2);
T_Energy[is]->Fit(myMIPs,"R");

EMCEnergy[is]->cd(16);
T_Energy[is]->Draw();
gPad->SetLogy();
EMCEnergy[is]->Print(EPSFile);
}
TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./Energy.exe ./Raw2ROOT/<filename>\n";
	cout<<"default output file name is <inputfilename.eps>"<<endl;	
	return;
}
void phrase_command(int argc,char *argv[]){
	if (argc<2){ helpinfo();
	     exit(0);
	}else {filename=(string) argv[1]; cout<<"START\n~~~~~*****~~~~~\nInput File : "<<filename<<endl;}
}

int main(int argc,char *argv[])
 {

phrase_command(argc,argv);
string Tt="CTest";
//string RootFileName=filename;
const char *TtreeName=(char *)(Tt.data());
const Text_t* RawRootFile=(Text_t*)(filename.data());
char *RawRootFile_c=(char *)(filename.data());
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(RawRootFile);
if (!f) 
  {
  f = new TFile(RawRootFile_c);
  }
 TTree *tree = (TTree*)gDirectory->Get(TtreeName);
 
 RootConvert xx;
 xx.RootConvert::Init(tree);
 xx.RootConvert::Loop(RawRootFile_c);
 cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
 return 1;
 }

