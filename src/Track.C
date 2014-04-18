#define Track_cxx
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
if(MaxADC[0]>150&&(MaxADC[1]>50||MaxADC[12]>50)&&MaxADC[13]>150)
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
string namesuffix="_Track.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[50]="./Track/";
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
TH1F *EMC_TraRes[Nlayer];
for (int nl=0;nl<Nlayer;nl++)
{
char cc[30];
sprintf(cc,"Layer%d_TrackResolution",nl);
EMC_TraRes[nl]=new TH1F(cc,cc,45,-1.5,1.5);
EMC_TraRes[nl]->SetXTitle("Res(2.7cm)");
EMC_TraRes[nl]->SetYTitle("Counts");
}
Double_t TraResPar[Nlayer][3];
for(int ig=0;ig<Nlayer;ig++)
{
memset(TraResPar[ig],0,sizeof(TraResPar[ig]));
}
gStyle->SetOptStat(111);
gStyle->SetOptFit(1111);
TH2F *TrackShow[24][2];
Int_t nShow[2]={0,0};
Double_t Par[2];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//LOOP START
cout<<"Loop root file ..."<<endl;
int BadTrack=0;
int Tevent=0;
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
//cout<<"222222222222222222222222"<<endl;
//track for hit position
TF1 *linear=new TF1("linear","[0]*x+[1]",-1,23);
TH2F *mytrack[2];
//mytrack[0]=new TH2F("mytrackX","matrackX",160,-1,15,240,-1,23);
//mytrack[1]=new TH2F("mytrackY","mytrackY",160,-1,15,240,-1,23);
mytrack[0]=new TH2F("mytrackX","matrackX",240,-1,23,160,-1,15);
mytrack[1]=new TH2F("mytrackY","mytrackY",240,-1,23,160,-1,15);
for(int il=0;il<Nlayer;il++)
{//cout<<"111111111111111111111111"<<endl;
//mytrack->Fill(FitTrack[il],il);
if(il%2==0)
//mytrack[1]->Fill(il+1,FitTrack[il]);//will be helpful for case of vertical tracks
mytrack[1]->Fill(FitTrack[il],il+1);//will be helpful for case of vertical tracks
else
//mytrack[0]->Fill(il+1,FitTrack[il]);
mytrack[0]->Fill(FitTrack[il],il+1);
//cout<<"222222222222222222222222"<<endl;
}
Par[0]=0;
Par[1]=0;
for(int idim=0;idim<2;idim++){
mytrack[idim]->Fit(linear,"RQ0");
linear->GetParameters(Par);
linear->SetParameters(Par);
mytrack[idim]->Fit(linear,"RQ0");
linear->GetParameters(Par);
Double_t ChiS=linear->GetChisquare();
 if(ChiS<10000){
    for(int il=0;il<Nlayer;il++){
      if(il%2!=idim){
      Double_t Res=0;
      //Res=(il+1)*Par[0]+Par[1]-FitTrack[il];
      Res=((il+1)-Par[1])/Par[0]-FitTrack[il];
      EMC_TraRes[il]->Fill(Res);
    }
//    cout<<FitTrack[il]<<" "<<endl;
    }
if(nShow[idim]<24)
{
char cS[30];
if(idim==0)
sprintf(cS,"Event%d_Track_ZX",Tevent);
if(idim==1)
sprintf(cS,"Event%d_Track_ZY",Tevent);
TrackShow[nShow[idim]][idim]=(TH2F*)mytrack[idim]->Clone();
TrackShow[nShow[idim]][idim]->SetName(cS);
TrackShow[nShow[idim]][idim]->SetTitle(cS);
nShow[idim]++;
}}
}
delete mytrack[0];
delete mytrack[1];
delete linear;
}
else 
{BadTrack++;
}
Tevent++;
}
//Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Track Fit and Draw
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());


TF1 *linear=new TF1("linear","[0]*x+[1]",-1,24);
TCanvas *trac[2];
for(int id=0;id<2;id++){
char ss[30];
sprintf(ss,"dimension%d",id);
trac[id]=new TCanvas(ss,ss,0,0,600,800);
trac[id]->Divide(4,6);
trac[id]->SetGrid();
for(int ishow=0;ishow<nShow[id];ishow++)
{
trac[id]->cd(ishow+1);
TrackShow[ishow][id]->Draw();
//TrackShow[ishow][id]->SetGridX();
//TrackShow[ishow][id]->SetGridY();
TrackShow[ishow][id]->SetMarkerStyle(25);
TrackShow[ishow][id]->SetMarkerSize(0.2);
TrackShow[ishow][id]->SetMarkerColor(kBlue);
TrackShow[ishow][id]->Fit(linear,"RQ");
}
trac[id]->Print(EPSFile);
}


TF1 *mygaus=new TF1("mygaus","gaus",-2,2);
TCanvas *B=new TCanvas("B","B",0,0,600,800);
B->Divide(4,4);
for(int il=0;il<Nlayer;il++)
{
B->cd(il+1);
EMC_TraRes[il]->Draw();
EMC_TraRes[il]->Fit(mygaus,"QR");
}
B->Print(EPSFile);


TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*ofstream outputpar;
outputpar.open("./TrackV2/TrackV2.txt");
if(!outputpar.good())
{
cout<<"can not open out put fitted par file :./Track/Track_LangauPar.txt"<<endl;
exit(1);
}
for(int i=0;i<nGID;i++)
{
if(i%3==2)
{
//const int nLangaus=4;
//Double_t TrackPar[nGID][nLangaus];
for(int j=0;j<nLangaus;j++)
{
outputpar<<TrackPar[i][j]<<" ";
}
outputpar<<"\n";
}
}
outputpar.close();
*/
cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;
}//Loop end 




//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./Track.exe ./Raw2ROOT/<filename>\n";
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

