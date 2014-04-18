#define DyRel_cxx
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
using namespace std;
using namespace cosmictest;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//BGO constants information
//Plane(1-7), Layer(1-14), Bar(1-44) ,Side(0,1),Dimension(0,1),FEEcard Number(1-16),FEEchannel Number(144)
//const int Nlayer=7;
//const int Nbar=24;
//const int Ndy=3;
//const int Nside=1;
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
if((MaxADC[0]>200||MaxADC[1]>200)||(MaxADC[4]>200||MaxADC[5]>200))
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
string namesuffix="_DyRel.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[80]="./DyRel/";
char *EPSFile=strcat(FilePath,outputfile.c_str());
cout<<"Output file: "<<EPSFile<<endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read pedestal data
cout<<"Read Pedestal..."<<endl;

const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
const int nGID=(Nlayer*Nside)*Nbar*Ndy;
Double_t PedestalPar[nGID][nDyPar];
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
TH2F *EMC_DyRel[nGID];
for (int ng=0;ng<nGID;ng++)
{
       // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
Int_t idy=ng%Ndy;//0:dy8/dy2 1:dy5/dy2 2:dy8/dy5
Int_t ibar=((int)(ng/Ndy))%Nbar;
Int_t iside=((int)(ng/Nbar/Ndy))%Nside;
Int_t ilayer=(int)(ng/Nside/Nbar/Ndy);
char cc[30];
if(idy==0)
sprintf(cc,"Layer%d_Side%d_Bar%d_Dy8/dy2",ilayer+1,iside,ibar);
else if(idy==1)
sprintf(cc,"Layer%d_Side%d_Bar%d_Dy5/dy2",ilayer+1,iside,ibar);
else if(idy==2)
sprintf(cc,"Layer%d_Side%d_Bar%d_Dy8/dy5",ilayer+1,iside,ibar);
EMC_DyRel[ng]=new TH2F(cc,cc,500,0,1000,3000,0,16500);
}
 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//LOOP START
cout<<"Loop root file ..."<<endl;
int BadTrack=0;
int Tevent=0;
Double_t FitRange[nGID][2];
for(int ig=0;ig<nGID;ig++)
{
memset(FitRange[ig],0,sizeof(FitRange[ig]));
}

if (fChain == 0) {cout<<"fChain=0,return!";return;}
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) {cout<<"ientry="<<ientry<<",break";break;}
  nb = fChain->GetEntry(jentry);   
  if(jentry%10000==0)
  cout<<"Total Events:"<<jentry<<endl;

 Int_t FitTrack[15];//
//void TrackSeek(Int_t nHits,Int_t *Layer,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][2],Int_t *Dy,Int_t *Side,Int_t nGID,Double_t FitTrack[15],Double_t FitADC[14])
 TrackSeek(nHits,Layer,Bar,ADC,PedestalPar,Dy,Side,GID,FitTrack);
if(FitTrack[14]==1)
{
Double_t dynode[nGID];
memset(dynode,0,sizeof(dynode));
for(int i=0;i<nHits;i++)
{
    if (ADC[i]>(PedestalPar[GID[i]][0]+2.5*PedestalPar[GID[i]][1]))
      dynode[GID[i]]=ADC[i]-PedestalPar[GID[i]][0];
}
for(int i=0;i<nGID;i++)
if(i%3==2)
{   EMC_DyRel[i]->Fill(dynode[i-1],dynode[i]);
    EMC_DyRel[i-1]->Fill(dynode[i-2],dynode[i-1]);
    
    if(dynode[i]>3000&&dynode[i]<4500&&FitRange[i][0]<50)    
    FitRange[i][0]=dynode[i-1];
    else if(dynode[i]>11000&&dynode[i]<12000&&FitRange[i][1]<50)    
    FitRange[i][1]=dynode[i-1];
}
}
else 
{BadTrack++;
}
Tevent++;
}
//Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DyRel Fit and Draw
TF1 *myDyRel=new TF1("myDyRel","[0]+[1]*x",0,2500);
myDyRel->SetParNames("Intercept","Slope");
const int nLinear=2;
Double_t DyRelPar[nGID][nLinear];
for(int ig=0;ig<nGID;ig++)
{
memset(DyRelPar[ig],0,sizeof(DyRelPar[ig]));
}
gStyle->SetOptStat(111);
gStyle->SetOptFit(1111);

//
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

TCanvas *EMCDyRel[Nlayer][Nside];
for (int il=0;il<Nlayer;il++)
for(int is=0;is<Nside;is++)
{
char d[50];
sprintf(d,"Layer%d_Side%d_Dy8/Dy5",il,is);
EMCDyRel[il][is]=new TCanvas(d,d,0,0,600,800);
EMCDyRel[il][is]->Divide(4,6);

for(int ib=0;ib<Nbar;ib++)
{
Int_t iGID=((il*Nside+is)*Nbar+ib)*Ndy+2;//idy=2;
EMCDyRel[il][is]->cd(ib+1);
myDyRel->SetRange(FitRange[iGID][0],FitRange[iGID][1]);
EMC_DyRel[iGID]->Draw();
myDyRel->SetLineColor(2);
if (FitRange[iGID][0]*FitRange[iGID][1]!=0)
  {
  EMC_DyRel[iGID]->Fit(myDyRel,"R");
  myDyRel->GetParameters(DyRelPar[iGID]);
  }
  else 
  memset(DyRelPar[iGID],0,sizeof(DyRelPar[iGID]));
}
EMCDyRel[il][is]->Print(EPSFile);
}

TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ofstream outputpar;
outputpar.open("./DyRel/DyRel.txt");
if(!outputpar.good())
{
cout<<"can not open out put fitted par file :./DyRel/DyRel.txt"<<endl;
exit(1);
}
for(int i=0;i<nGID;i++)
{
if(i%3==2)
{
//const int nLangaus=4;
//Double_t DyRelPar[nGID][nLangaus];
for(int j=0;j<2;j++)
{
outputpar<<DyRelPar[i][j]<<" ";
}
outputpar<<"\n";
}
}
outputpar.close();
cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./DyRel.exe ./Raw2ROOT/<filename>\n";
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

//string RootFileName=filename;
phrase_command(argc,argv);
string Tt="CTest";

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

