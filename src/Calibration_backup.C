#define Calibration_cxx
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
#include <RootConvert.h>
#include "TFile.h"
#include <TF1.h>
#include "TInterpreter.h"
#include <TStyle.h>
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"



using namespace std;
using namespace cosmictest;

void Order(Float_t d[],Float_t *source){
//////////////
Float_t t=0;
Int_t   x=0;
for(int i=0;i<13;i++)
d[i]=source[i];

for(int i=0;i<13;i++){
t=d[i];x=i;
for(int j=i+1;j<13;j++){
if(t>=d[j])
{t=d[j];x=j;
}
}
d[x]=d[i];d[i]=t;
}
for(int i=0;i<13;i++)
cout<<d[i]<<"  ";

cout<<endl;

}

void RootConvert::Loop(char *RootFileName)
{

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//output .root .eps files
//string namesuffix="_pedestal.eps";
string namesuffix="_calibration";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);

int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
//pedestal par cout 
string TXTFile0="./Calibration/"+out+".txt";
outputfile=out+".eps";

  //char FilePath[80]="./Calibration/";
  
  string  EPSFile="./Calibration/"+outputfile;
// string TXTFile=(string)EPSFile-".eps"+".txt";

string TXTFile1="./Calibration/Calibration.txt";
cout<<"Output file 1 "<<EPSFile<<endl;
cout<<"Output file 2 "<<TXTFile0<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define the Calibration output
const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
Double_t CalibrationPar[nGID][nDyPar];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define the TH1F for pedestal 
Float_t Cal_DACcode[15]={160,320,480,640,800,960,1120,1280,1440,1600,1760,1920,2080,2240,2400};
Float_t Cal_Vol[15];
for(int i=0;i<15;i++){
Cal_Vol[i]=0.75651*Cal_DACcode[i]+2.64069;
}
TH1F *EMC_Cal[nGID];
TH2F *FEE_G[nGID];
for (int ng=0;ng<nGID;ng++)
{
       // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
Int_t idy=ng%Ndy;
Int_t ibar=((int)(ng/Ndy))%Nbar;
Int_t iside=((int)(ng/Nbar/Ndy))%Nside;
Int_t ilayer=(int)(ng/Nside/Nbar/Ndy);
char cc[40];
sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_cal",ilayer,iside,ibar,(idy*3+2));
char dd[40];
sprintf(dd,"Layer%d_Side%d_Bar%d_Dy%d_FEE",ilayer,iside,ibar,(idy*3+2));
char ee[40];
sprintf(dd,"Layer%d_Side%d_Bar%d_Dy%d_FEE;mV;ADC channels",ilayer,iside,ibar,(idy*3+2));
EMC_Cal[ng]=new TH1F(cc,cc,8500,-1000,16000);
FEE_G[ng]=new TH2F(dd,ee,250,0,2500,1700,-1000,16000);
}

TF1 *mylinear=new TF1("mylinear","[0]*x+[1]",100,800);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Loop start
cout<<"Loop Root File..."<<endl;
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
cout<<"~~~~~~~~~~~~~~~"<<nentries<<endl;
Long64_t nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   
 
//if(time<19370)
//continue;
//else if (time>24270)
//break;
  if(jentry%10000==0)
  cout<<"Total Events:"<<jentry<<endl;
  for(int i=0; i<nHits;i++)
  EMC_Cal[GID[i]]->Fill(ADC[i]);
 }//single event Loop end


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Calibration Fit and Draw
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

Double_t par[2];
TSpectrum *s=new TSpectrum(15,50);
Int_t nfound =0;
//printf("Found %d candidate peaks \n",nfound);
Float_t *xpeaks;

TCanvas *EMC_Ped[Nlayer][Nside][Ndy];
for(int ilayer=0;ilayer<Nlayer;ilayer++)
for(int iside=0;iside<Nside;iside++)
{
char d[3][50];
//Dy2,Dy5,Dy8
for(int idy=0;idy<Ndy;idy++)
{
sprintf(d[idy],"Layer%d_Side%d_Dy%d",ilayer+1,iside,idy*3+2);
EMC_Ped[ilayer][iside][idy]=new TCanvas(d[idy],d[idy],0,0,600,800);
EMC_Ped[ilayer][iside][idy]->Divide(4,6);

for(int ibar=0;ibar<Nbar;ibar++)
{
Int_t iGID=((ilayer*Nside+iside)*Nbar+ibar)*Ndy+idy;
EMC_Ped[ilayer][iside][idy]->cd(ibar+1);
nfound=s->Search(EMC_Cal[iGID],22,"",0.10);
xpeaks=s->GetPositionX();
Float_t pp[15];
Order(pp,xpeaks);

for(int npeak=0;npeak<13;npeak++)
FEE_G[iGID]->Fill(Cal_Vol[npeak],pp[npeak]);

mylinear->SetRange(50,1500);
FEE_G[iGID]->Draw();
FEE_G[iGID]->SetMarkerStyle(24);
FEE_G[iGID]->SetMarkerSize(0.5);
FEE_G[iGID]->Fit(mylinear,"R");
mylinear->GetParameters(par);
//extract the Calibration parameters
CalibrationPar[iGID][0]=par[0];
CalibrationPar[iGID][1]=par[1];
}
EMC_Ped[ilayer][iside][idy]->Print(EPSFile.c_str());
}

}
TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//cout CalibrationPar
ofstream PedPar0;
PedPar0.open(TXTFile0.c_str());
if (!PedPar0.good())
{cout<<"Can not open Calibration TXTFlie Name File!!!"<<endl;
exit(0);
}
//Double_t CalibrationPar[nTBar][nDyPar];
for(int ig=0;ig<nGID;ig++)
{
PedPar0<<CalibrationPar[ig][0]<<"  "<<CalibrationPar[ig][1]<<" ";
PedPar0<<"\n";
}
PedPar0.close();

//cout CalibrationPar
ofstream PedPar1;
PedPar1.open(TXTFile1.c_str());
if (!PedPar1.good())
{cout<<"Can not open Calibration TXTFlie Name File!!!"<<endl;
exit(0);
}
//Double_t CalibrationPar[nTBar][nDyPar];
for(int ig=0;ig<nGID;ig++)
{
PedPar1<<CalibrationPar[ig][0]<<"  "<<CalibrationPar[ig][1]<<" ";
PedPar1<<"\n";
}
PedPar1.close();

//filename
ofstream name;

name.open("./Calibration/filename");
if(!name.good())
{cout<<"Can not open namefile !"<<endl;exit(1);}
name<<RootFileName;
name.close();

cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./Calibration.exe ./Raw2ROOT/<filename>\n";
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
 return 1;

 }

