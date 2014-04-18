#define Att_cxx
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
string namesuffix="_AttV2.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[50]="./AttV2/";
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
//TH2F *EMC_Att[nGbar];
TProfile *EMC_Att[nGbar];
for (int ng=0;ng<nGbar;ng++)
{
       // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
Int_t ibar=ng%Nbar;
Int_t ilayer=(int)(ng/Nbar);
char cc[30];
sprintf(cc,"Layer%d_Bar%d_Att",ilayer,ibar);
//EMC_Att[ng]=new TH2F(cc,cc,600,0,60,100,-5,5);
EMC_Att[ng]=new TProfile(cc,cc,22,0,60);
}
 
/*TF1 *myAtt=new TF1("myAtt",langaufun,0,2500,4);
myAtt->SetParNames("Width","MP","Area","GSigma");
const int nLangaus=4;
Double_t AttPar[nGID][nLangaus];
for(int ig=0;ig<nGID;ig++)
{
memset(AttPar[ig],0,sizeof(AttPar[ig]));
}*/
gStyle->SetOptStat(111);
gStyle->SetOptFit(0111);
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
if(FitTrack[14]==1)
{
//track for hit position
TF1 *linear=new TF1("linear","[0]*x+[1]",-1,24);
TH2F *mytrack[2];
mytrack[0]=new TH2F("mytrack_X","mytrack_X",160,-1,15,240,-1,23);
mytrack[1]=new TH2F("mytrack_Y","mytrack_Y",160,-1,15,240,-1,23);
for(int il=0;il<Nlayer;il++)
{
//mytrack->Fill(FitTrack[il],il);
if(il%2==0)
mytrack[1]->Fill(il+1,FitTrack[il]);//will be helpful for case of vertical tracks
else
mytrack[0]->Fill(il+1,FitTrack[il]);
}
Par[0]=0;Par[1]=0;
Double_t HitPos[Nlayer];//Hitted position of BGO
memset(HitPos,0,sizeof(HitPos));
for(int idim=0;idim<2;idim++){
mytrack[idim]->Fit(linear,"RQ0");
linear->GetParameters(Par);
linear->SetParameters(Par);
mytrack[idim]->Fit(linear,"RQ0");
linear->GetParameters(Par);
Double_t ChiS=linear->GetChisquare();
  if(ChiS<10){
    for(int il=0;il<Nlayer;il++){
      if(il%2==idim)
      HitPos[il]=((Par[0]*(il+1)+Par[1])-0.5)*2.75;//unit :cm
//cout<<"HitPosition:"<<HitPos[il]<<endl;
    }
  }
  else memset(HitPos,0,sizeof(HitPos));
}
delete linear;
delete mytrack[0];
delete mytrack[1];
//log(A0/A1)
Double_t Am[Nlayer][Nside]={{'0'},{'0'}};
for(int i=0;i<nHits;i++){
if(FitTrack[Layer[i]]==Bar[i]&&Dy[i]==8)//both sides
    if (ADC[i]>(PedestalPar[GID[i]][0]+3*PedestalPar[GID[i]][1]))
    Am[Layer[i]][Side[i]]=ADC[i]-PedestalPar[GID[i]][0];
}
//Fill EMC_Att
for(int il=0;il<Nlayer;il++){
int ib=FitTrack[il];
int iGbar=il*Nbar+ib;
//cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//cout<<"Am_side0:"<<Am[il][0]<<" "<<MIPsPar[(il*Nside+0)*Nbar+ib][1]<<endl;
//cout<<"Am_side1:"<<Am[il][1]<<" "<<MIPsPar[(il*Nside+1)*Nbar+ib][1]<<endl;
Am[il][0]=Am[il][0]/MIPsPar[(il*Nside+0)*Nbar+ib][1];
Am[il][1]=Am[il][1]/MIPsPar[(il*Nside+1)*Nbar+ib][1];

//cout<<"HitPosition:"<<HitPos[il]<<endl;
if(Am[il][0]!=0&&Am[il][1]!=0&&HitPos[il]!=0)
{
//cout<<"HitPosition:"<<HitPos[il]<<endl;
EMC_Att[iGbar]->Fill(HitPos[il],TMath::Log(Am[il][0]/Am[il][1]));
}
}
}
else 
{BadTrack++;
}
Tevent++;
}
//Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Att Fit and Draw
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

TCanvas *EMCAtt[Nlayer];
for(int il=0;il<Nlayer;il++)
{
char d[50];
sprintf(d,"Layer%d_AttCff",il);
EMCAtt[il]=new TCanvas(d,d,0,0,600,800);
EMCAtt[il]->Divide(4,6);
for(int ib=0;ib<Nbar;ib++)
{
Int_t iGbar=il*Nbar+ib;
EMCAtt[il]->cd(ib+1);
EMC_Att[iGbar]->Draw();
}
EMCAtt[il]->Print(EPSFile);
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
	cout<<"Usage is ./Att.exe ./Raw2ROOT/<filename>\n";
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

