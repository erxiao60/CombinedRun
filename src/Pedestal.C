#define Pedestal_cxx
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

using namespace std;
using namespace cosmictest;
void RootConvert::Loop(char *RootFileName)
{

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//output .root .eps files
//string namesuffix="_pedestal.eps";
string namesuffix="_pedestal";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);

int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
//pedestal par cout 
string TXTFile0="./Pedestal/"+out+".txt";
outputfile=out+".eps";

  char FilePath[80]="./Pedestal/";
  char *EPSFile=strcat(FilePath,outputfile.c_str());
// string TXTFile=(string)EPSFile-".eps"+".txt";

string TXTFile1="./Pedestal/Pedestal.txt";
cout<<"Output file 1 "<<EPSFile<<endl;
cout<<"Output file 2 "<<TXTFile0<<endl;
//  TFile *outFile=new TFile(ROOTFile,"RECREATE");
//  if (outFile==(TFile*) NULL)
//  {
//    std::cerr<<"Error:Could not Open Conversion ROOT File!"<<endl;
//    exit(1);
//  }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//BGO constants information
//Plane(1-7), Layer(1-14), Bar(1-44) ,Side(0,1),Dimension(0,1),FEEcard Number(1-16),FEEchannel Number(144)
/*
ifstream BGOinfor;
BGOinfor.open("./map/BGO.infor");
if(!BGOinfor.good())
{
cout<<"~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Can not open BGO constant information file!"<<endl;
exit(1);
}
int ConBGO[5];
char infnote[2][80];
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"~~~~BGOECAL Cosmic Test~~~~"<<endl;
BGOinfor.getline(infnote[0],80);
BGOinfor.getline(infnote[1],80);
cout<<infnote[1]<<endl;
for(int inf=0;inf<5;inf++)
{
BGOinfor>>ConBGO[inf];
cout<<ConBGO[inf]<<"  ";
}
cout<<"\n";
const int Nlayer=ConBGO[0];
const int Nbar=ConBGO[1];
const int Ndy=ConBGO[2];
const int Nside=ConBGO[3];
//const int N_FEE=ConBGO[4];
BGOinfor.close();
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define the Pedestal output
const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
//const int nGID=(Nlayer*Nside)*Nbar*Ndy;
Double_t PedestalPar[nGID][nDyPar];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define the TH1F for pedestal 
TH1F *EMC_Cal[nGID];
for (int ng=0;ng<nGID;ng++)
{
       // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
Int_t idy=ng%Ndy;
Int_t ibar=((int)(ng/Ndy))%Nbar;
Int_t iside=((int)(ng/Nbar/Ndy))%Nside;
Int_t ilayer=(int)(ng/Nside/Nbar/Ndy);
char cc[30];
//for quadrant 1 ,3
sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_Ped",ilayer,iside,ibar,(idy*3+2));
EMC_Cal[ng]=new TH1F(cc,cc,800,-400,400);
}

TF1 *PedFit=new TF1("PedFit","gaus",100,800);

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
//    cout<<"~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//    cout<<"Tevent:"<<jentry<<" nHits:"<<nHits<<endl;
  for(int i=0; i<nHits;i++){
  EMC_Cal[GID[i]]->Fill(ADC[i]);
//  if(Dy[i]==8){//both sides
//    if (ADC[i]==0)
//    cout<<"Layer:"<<Layer[i]<<" Bar:"<<Bar[i]<<" Side:"<<Side[i]<<" ADC:"<<ADC[i]<<endl;
//}
}
}
gStyle->SetOptFit(1111);
gStyle->SetOptStat(111);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Pedestal Fit and Draw
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

Double_t par[3];
TSpectrum *s=new TSpectrum(1,30);
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
PedFit->SetRange(xpeaks[0]-20,xpeaks[0]+20);
PedFit->SetParameters(1000,xpeaks[0]);
EMC_Cal[iGID]->Fit(PedFit,"R0");
PedFit->GetParameters(par);
PedFit->SetRange(par[1]-3*par[2],par[1]+3*par[2]);
PedFit->SetParameters(par);
EMC_Cal[iGID]->SetAxisRange(par[1]-100,par[1]+100,"X");

EMC_Cal[iGID]->Fit(PedFit,"R0");
PedFit->GetParameters(par);
PedFit->SetRange(par[1]-2.5*par[2],par[1]+2.5*par[2]);
PedFit->SetParameters(par);
PedFit->SetLineColor(2);
EMC_Cal[iGID]->Draw();
EMC_Cal[iGID]->Fit(PedFit,"R");

//extract the Pedestal parameters
PedFit->GetParameters(par);
PedestalPar[iGID][0]=par[1];
PedestalPar[iGID][1]=par[2];
}
EMC_Ped[ilayer][iside][idy]->Print(EPSFile);
}

}
TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//cout PedestalPar
ofstream PedPar0;
PedPar0.open(TXTFile0.c_str());
if (!PedPar0.good())
{cout<<"Can not open Pedestal TXTFlie Name File!!!"<<endl;
exit(0);
}
//Double_t PedestalPar[nTBar][nDyPar];
for(int ig=0;ig<nGID;ig++)
{
PedPar0<<PedestalPar[ig][0]<<"  "<<PedestalPar[ig][1]<<" ";
PedPar0<<"\n";
}
PedPar0.close();

//cout PedestalPar
ofstream PedPar1;
PedPar1.open(TXTFile1.c_str());
if (!PedPar1.good())
{cout<<"Can not open Pedestal TXTFlie Name File!!!"<<endl;
exit(0);
}
//Double_t PedestalPar[nTBar][nDyPar];
for(int ig=0;ig<nGID;ig++)
{
PedPar1<<PedestalPar[ig][0]<<"  "<<PedestalPar[ig][1]<<" ";
PedPar1<<"\n";
}
PedPar1.close();

//filename
ofstream name;

name.open("./Pedestal/filename");
if(!name.good())
{cout<<"Can not open namefile !"<<endl;exit(1);}
name<<RootFileName;
name.close();

cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./Pedestal.exe ./Raw2ROOT/<filename>\n";
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

