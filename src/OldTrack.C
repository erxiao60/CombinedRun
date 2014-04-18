#define Track_cxx
#include "Track.h"
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
TF1 *Tracks=new TF1("Tracks","[0]*x+[1]",0,15); 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//BGO constants information
//Plane(1-7), Layer(1-14), Bar(1-44) ,Side(0,1),Dimension(0,1),FEEcard Number(1-16),FEEchannel Number(144)
const int Nplane=4;
const int Nlayer=8;
const int Nbar=44;
const int Ndy=3;
const int Nside=2;
const int Ndimension=2;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TrackSeek(Double_t EneCor[][Nbar/2],Double_t SecCoo[][Nbar/2],Double_t *FirCooY,Double_t *SecCooX,int *nfalseL,Int_t *MaxBar)
{
  Double_t pitch=37/12;//cm
  Double_t MaxADC[Nlayer];
//Int_t   MaxBar[Nlayer];
//  Double_t FirCooY[Nlayer];
//  Double_t SecCooX[Nlayer];
  memset(MaxADC,0,sizeof(MaxADC));
  memset(MaxBar,0,sizeof(MaxBar));
  memset(FirCooY,0,sizeof(FirCooY));
  memset(SecCooX,0,sizeof(SecCooX));
 for(int il=0;il<Nlayer;il++)
 for(int ib=0;ib<Nbar/2;ib++)
 if(EneCor[il][ib]>MaxADC[il])
 {
 MaxADC[il]=EneCor[il][ib];
 MaxBar[il]=ib;
 } 
 
cout<<"MaxADC:"<<endl; 
 for(int il=0;il<Nlayer;il++)
{
cout<<MaxADC[il]<<" ";
if(MaxADC[il]<500)
nfalseL[0]++;
}
 for(int il=0;il<Nlayer;il++)
 {
 Int_t b=MaxBar[il];
 Double_t TE=EneCor[il][b];
 Double_t TEB=TE*b;
 if(MaxBar[il]!=0)
{ TE+=EneCor[il][b-1];
  TEB+=EneCor[il][b-1]*(b-1);
}
 if(MaxBar[il]<Nbar/2-1)
{ TE+=EneCor[il][b+1];
  TEB+=EneCor[il][b+1]*(b+1);
}
FirCooY[il]=(TEB/TE+0.5)*pitch;
SecCooX[il]=SecCoo[il][b];
}
}
void Track::Loop(char *RootFileName)
{

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
const int nTBar=Nlayer*Nbar;//12 layers * 13bars
const int nDyPar=Ndy*2;//Dy2,5,8 mean and sigma (gaus fitted)
Double_t PedestalPar[nTBar][6];
ifstream PedPar;
PedPar.open("./Pedestal/Pedestal.txt");
if (!PedPar.good())
{cout<<"Can not open Pedestal TXTFlie  File!!!"<<endl;
exit(0);
}
for(int iTBar=0;iTBar<nTBar;iTBar++)
for(int iDyPar=0;iDyPar<nDyPar;iDyPar++)
PedPar>>PedestalPar[iTBar][iDyPar];
PedPar.close();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read Attenuate data
Double_t AttFitPar[Nlayer][Nbar][2];
ifstream attpar;
attpar.open("./Att/Att.txt");
if(!attpar.good())
{
cout<<"can not open out put fitted par file :./Att/Att.txt"<<endl;
exit(1);
}
for(int i=0;i<Nlayer;i++)
for(int j=0;j<Nbar;j++)
{
for(int n=0;n<2;n++)
attpar>>AttFitPar[i][j][n];
if(AttFitPar[i][j][1]>200||AttFitPar[i][j][1]<100)
AttFitPar[i][j][1]=140;
}
attpar.close();

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//std::vector<std::vector<TH2D*> > CosTrack;
std::vector<TH1D*>  SpaResFir;
std::vector<TH1D*>  SpaResSec;
char cc1[Nlayer][30];
char cc2[Nlayer][30];
for(Int_t ilayer=0;ilayer<Nlayer;ilayer++)
 { 
 sprintf(cc1[ilayer],"Layer%d_SpaceResFir",(unsigned)(ilayer+1));
 SpaResFir.push_back(new TH1D(cc1[ilayer],cc1[ilayer],200,-5,5));
 sprintf(cc2[ilayer],"Layer%d_SpaceResSec",(unsigned)(ilayer+1));
 SpaResSec.push_back(new TH1D(cc2[ilayer],cc2[ilayer],200,-5,5));
 }

//TH2D *Coeff[Nlayer][Nbar];
TProfile *Coeff[Nlayer][Nbar];
for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar/2;ib++)//Bar 0===Bar 44 ,1==43 ,A,B end of a same BGO bar
{
char coe[2][50];
sprintf(coe[0],"Layer%d_BGO%d_Side%d_Coe",il+1,ib+1,0);
sprintf(coe[1],"Layer%d_BGO%d_Side%d_Coe",il+1,ib+1,1);
//Coeff[il][ib]=new TH2D(coe[0],coe[0],120,0,60,250,300,5000);
//Coeff[il][43-ib]=new TH2D(coe[1],coe[1],120,0,60,250,300,5000);
Coeff[il][ib]=new TProfile(coe[0],coe[0],60,0,60);
Coeff[il][43-ib]=new TProfile(coe[1],coe[1],60,0,60);
}
TF1 *myfit=new TF1("myfit","[0]+[1]*x",0,60);
gStyle->SetOptStat(111);
gStyle->SetOptFit(0111);

TH2D *TrackShow[100];
TH2D *TrackShowFir[100];
for(int ts=0;ts<100;ts++)
{
char tt[30];
char tl[30];
sprintf(tt,"Track_%d",ts+1);
sprintf(tl,"TrackLow_%d",ts+1);
TrackShow[ts]=new TH2D(tt,"Track;x or y;z",200,-5,75,100,0,35);
TrackShowFir[ts]=new TH2D(tl,"Track;x or y;z",200,-5,75,100,0,35);
}

int nShow=0;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Loop start
int Tevent=0;
int BadTrack=0;
Double_t EneCor[Nlayer][Nbar/2];//Acorr :corrected Energy 
// Initialize 
for(int ii=0;ii<Nlayer;ii++)
for(int jj=0;jj<Nbar/2;jj++)
EneCor[ii][jj]=0;

Double_t SecCoo[Nlayer][Nbar/2];//Second Coordinate
// Initialize 
for(int ii=0;ii<Nlayer;ii++)
for(int jj=0;jj<Nbar/2;jj++)
SecCoo[ii][jj]=0;

Double_t Dynode[Nlayer][Nbar][Ndy];
for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar;ib++)  
for(int ir=0;ir<Ndy;ir++)
Dynode[il][ib][ir]=0;

if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   
  if(jentry%10000==0)
  cout<<"Total Events:"<<jentry<<endl;
//if(jentry==1000000)
//break;
  Tevent++;
 for(int i=0;i<nHits;i++)
      {
      int nTlayers=Plane[i]*2+Dimension[i];
      int nTbars=nTlayers*44+Bar[i]-1;
      int flag=(PMTdy[i]-2)/3;
//    cout<<"ADC:"<<ADC[i]<<" Plane:"<<Plane[i]<<" Dimension:"<<Dimension[i]<<" Bar:"<<Bar[i]<<endl;
      if (ADC[i]>(PedestalPar[nTbars][flag*2]+3*PedestalPar[nTbars][flag*2+1]))
        {
        
        Double_t ADCCutPed=ADC[i]-PedestalPar[nTbars][flag*2];
//      cout<<ADCCutPed<<endl;
   //     cout<<nTlayers<<" "<<nTbars<<" "<<flag<<endl;
        Dynode[nTlayers][Bar[i]-1][flag]=ADCCutPed;
        }
       }
  Double_t lambda;
  for(int il=0;il<Nlayer;il++)//using Dy8 only ,it is not precise!for preliminary results.
  for(int ib=0;ib<Nbar/2;ib++)  
  {
  lambda=(AttFitPar[il][ib][1]+AttFitPar[il][43-ib][1])/2;
  EneCor[il][ib]=sqrt(Dynode[il][ib][2]*Dynode[il][43-ib][2]*exp(60.0/lambda));
  SecCoo[il][ib]=(60-log(Dynode[il][ib][2]*AttFitPar[il][43-ib][0]/Dynode[il][43-ib][2]/AttFitPar[il][ib][0])*lambda)/2;
  }

  Double_t FirCooY[Nlayer];
  Double_t SecCooX[Nlayer];
  Int_t   MaxBar[Nlayer];
 int nfalseL=0;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"\n"<<"event:"<<jentry<<endl;
for(int l=0;l<Nlayer;l++)
cout<<FirCooY[l]<<" ";
cout<<"\n";
for(int l=0;l<Nlayer;l++)
cout<<SecCooX[l]<<" ";

 TrackSeek(EneCor,SecCoo,FirCooY,SecCooX,&nfalseL,MaxBar);
 if(nfalseL>0)
 {
 BadTrack++;
 continue;
 }
 int nPoint=0;
for(int ll=0;ll<Nlayer;ll++)
if(SecCooX[ll]<-5||SecCooX[ll]>65)
nPoint++;
if(nPoint>2)
 {
 BadTrack++;
 continue;
 }

//Fit the tracks


TH2D *Tr_x=new TH2D("Tr_x","Tr_x",300,-5,75,100,0,35);
TH2D *Tr_xFir=new TH2D("Tr_xFir","Tr_xFir",300,-5,75,100,0,35);
TH2D *Tr_y=new TH2D("Tr_y","Tr_y",300,-5,75,100,0,35);
TH2D *Tr_yFir=new TH2D("Tr_yFir","Tr_yFir",300,-5,75,100,0,35);
for(int jl=0;jl<Nlayer;jl++)
{
if(jl%2==0)
{Tr_x->Fill(SecCooX[jl],(jl+0.5)*3.5);
//Tr_y->Fill(FirCooY[jl],(jl+0.5)*3.5);
Tr_yFir->Fill(FirCooY[jl],(jl+0.5)*3.5);
}
else
{Tr_y->Fill(SecCooX[jl],(jl+0.5)*3.5);
//Tr_x->Fill(FirCooY[jl],(jl+0.5)*3.5);
Tr_xFir->Fill(FirCooY[jl],(jl+0.5)*3.5);
}
}
Double_t FitPar[2][2]={{0.,0.},{0.,0.}};
Double_t ChiSx=0.;
Double_t ChiSy=0.;
//Tr_x->Fit(myfit,"0");
Tr_xFir->Fit(myfit,"0");
myfit->GetParameters(FitPar[0]);
ChiSx=myfit->GetChisquare();
if(ChiSx<500&&nShow<100)
{
TrackShow[nShow]=(TH2D*)Tr_x->Clone();
TrackShowFir[nShow]=(TH2D*)Tr_xFir->Clone();
//TrackShow[nShow]->SetName("TrackShow");
nShow++;
}
//Tr_y->Fit(myfit,"0");
Tr_yFir->Fit(myfit,"0");
myfit->GetParameters(FitPar[1]);
ChiSy=myfit->GetChisquare();
delete Tr_x;
delete Tr_y;
delete Tr_xFir;//weighting
delete Tr_yFir;//weighting

//computing resolution
//Fill Coeff[][];
//std::vector<TH1D*>  SpaResFir;
//std::vector<TH1D*>  SpaResSec;
if(ChiSx<500&&ChiSy<500)
{
Double_t x,y,z;
x=0.;y=0.;z=0.;
for(int kl=0;kl<Nlayer;kl++)
{
 z=(kl+0.5)*3.5;
if(FitPar[0][1]!=0)
 x=(z-FitPar[0][0])/FitPar[0][1];
if(FitPar[1][1]!=0)
 y=(z-FitPar[1][0])/FitPar[1][1];
if(kl%2==0)
{
SpaResFir[kl]->Fill(y-FirCooY[kl]);
SpaResSec[kl]->Fill(x-SecCooX[kl]);

}
else
{
SpaResFir[kl]->Fill(x-FirCooY[kl]);
SpaResSec[kl]->Fill(y-SecCooX[kl]);
}
//Fill Coeff[][]
int hitbar=MaxBar[kl];
Double_t axis_x;
if(kl%2==0)
axis_x=x;
else
axis_x=y;
Coeff[kl][hitbar]->Fill(axis_x,Dynode[kl][hitbar][2]);
Coeff[kl][43-hitbar]->Fill(axis_x,Dynode[kl][43-hitbar][2]);
}
}
//reset members
for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar/2;ib++)  
EneCor[il][ib]=0.;

for(int ii=0;ii<Nlayer;ii++)
for(int jj=0;jj<Nbar/2;jj++)
SecCoo[ii][jj]=0;

for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar;ib++)  
for(int ir=0;ir<Ndy;ir++)
Dynode[il][ib][ir]=0;

}//single event Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Draw
TCanvas *SpaceFir=new TCanvas("SpaceFirst","SpaceFirst",0,0,600,800);
SpaceFir->Divide(3,3);
for(int l=0;l<Nlayer;l++)
{
SpaceFir->cd(l+1);
SpaResFir[l]->Draw();
}
SpaceFir->Print("./Track/resolution.eps(");
TCanvas *SpaceSec=new TCanvas("SpaceSecond","SpaceSecond",0,0,600,800);
SpaceSec->Divide(3,3);
for(int l=0;l<Nlayer;l++)
{
SpaceSec->cd(l+1);
SpaResSec[l]->Draw();
}
SpaceSec->Print("./Track/resolution.eps)");

TCanvas *TraShow=new TCanvas("S","S",0,0,600,800);
TraShow->Divide(4,6);
for(int ns=0;ns<24;ns++)
{
TraShow->cd(ns+1);
TrackShowFir[ns]->SetMarkerStyle(25);
TrackShowFir[ns]->SetMarkerSize(0.4);
TrackShowFir[ns]->SetMarkerColor(2);
TrackShowFir[ns]->Draw();
TrackShowFir[ns]->Fit(myfit);
TrackShow[ns]->SetMarkerStyle(23);
TrackShow[ns]->SetMarkerSize(0.4);
TrackShow[ns]->SetMarkerColor(4);
TrackShow[ns]->Draw("same");
}
TraShow->Print("./Track/TrackShow.eps");

TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
TCanvas *End=new TCanvas("End","End",0,0,600,800);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Attenuation coefficient
TH2D *CoeBar[Nside];//Bar 1=88 (bar+43*layer)
CoeBar[0]=new TH2D("CoeBar0","CoeVsBar_Side0;Tot_Bar;Coe(cm)",352,0,352,50,0,250);
CoeBar[1]=new TH2D("CoeBar1","CoeVsBar_Side1;Tot_Bar;Coe(cm)",352,0,352,50,0,250);
TH1D *CoeDis[Nside];
CoeDis[0]=new TH1D("CoeDis0","CoeDis_Side0;lambda(cm);counts",50,0,250);
CoeDis[1]=new TH1D("CoeDis1","CoeDis_Side1;lambda(cm);counts",50,0,250);
Double_t ChiSqr=0;
Double_t AttFunPar[Nlayer][Nbar][2];
TF1 *AttFun[2];
AttFun[0]=new TF1("AttFun0","[0]*exp(x/[1])",0,60);//37/12:BGO pitch
AttFun[1]=new TF1("AttFun1","[0]*exp((60-x)/[1])",0,60);//37/12:BGO pitch
Start->Print("./Att/Coefficientinteraction.eps[");
TCanvas *ccc[Nlayer][2];//2 layers, 2 sides
for(int il=0;il<Nlayer;il++)
for(int is=0;is<2;is++)
{char CoeC[30];
sprintf(CoeC,"Layer%d_Side%d",il+1,is+1);
ccc[il][is]=new TCanvas(CoeC,CoeC,0,0,600,800);
ccc[il][is]->Divide(4,6);
for(int ib=0;ib<Nbar/2;ib++)
{
ccc[il][is]->cd(ib+1);
int ttbar;
if(is==0) ttbar=ib;
if(is==1) ttbar=43-ib;
Coeff[il][ttbar]->Draw();
Coeff[il][ttbar]->SetMarkerStyle(23);
Coeff[il][ttbar]->SetMarkerSize(0.3);
Coeff[il][ttbar]->SetMarkerColor(4);
AttFun[is]->SetParameters(1200,-130);
AttFun[is]->SetLineColor(2);
Coeff[il][ttbar]->Fit(AttFun[is],"R");
ChiSqr=AttFun[is]->GetChisquare();
if(ChiSqr>5)
{AttFun[is]->SetRange(5,55);
AttFun[is]->SetLineColor(2);
Coeff[il][ttbar]->Fit(AttFun[is],"R");
}
AttFun[is]->GetParameters(AttFunPar[il][ttbar]);
AttFunPar[il][ttbar][1]=0-AttFunPar[il][ttbar][1];

CoeBar[is]->Fill(il*44+ttbar,AttFunPar[il][ttbar][1]);
CoeDis[is]->Fill(AttFunPar[il][ttbar][1]);
}
ccc[il][is]->Print("./Att/Coefficientinteraction.eps");
}
End->Print("./Att/Coefficientinteraction.eps]");
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~
//show coefficient distribution
TCanvas *A=new TCanvas("Coe_Dis","Coe_Dis",0,0,600,600);
A->Divide(2,3);
for(int is=0;is<Nside;is++)
{
A->cd(is*2+1);
CoeBar[is]->Draw();
CoeBar[is]->SetMarkerStyle(23);
CoeBar[is]->SetMarkerSize(0.3);
CoeBar[is]->SetMarkerColor(2*is+2);
A->cd(is*2+2);
CoeDis[is]->Draw();
}
gStyle->SetOptStat(0000);
A->cd(5);
CoeBar[0]->Draw();
CoeBar[0]->SetMarkerStyle(23);
CoeBar[0]->SetMarkerSize(0.3);
CoeBar[0]->SetMarkerColor(2);
CoeBar[1]->Draw("SAME");
CoeBar[1]->SetMarkerStyle(24);
CoeBar[1]->SetMarkerSize(0.3);
CoeBar[1]->SetMarkerColor(4);
A->cd(6);
CoeDis[0]->Draw();
CoeDis[0]->SetLineColor(2);
CoeDis[1]->Draw("SAME");
CoeDis[1]->SetLineColor(4);

A->Print("./Att/Dis_CoeIn.eps","Portrait");
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ofstream outputpar;
outputpar.open("./Att/AttInteraction.txt");
if(!outputpar.good())
{
cout<<"can not open out put fitted par file :./Att/AttInteraction.txt"<<endl;
exit(1);
}
for(int i=0;i<Nlayer;i++)
for(int j=0;j<Nbar/2;j++)
{
for(int n=0;n<2;n++)
outputpar<<AttFunPar[i][j][n]<<" ";
for(int n=0;n<2;n++)
outputpar<<AttFunPar[i][43-j][n]<<" ";
outputpar<<"\n";
}
outputpar.close();

cout<<"Total event:"<<Tevent<<"\n"<<"BadTrack:"<<BadTrack<<endl;
//Loop end 
}

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
string Tt="VA32Data";
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
 
 Track xx;
 xx.Track::Init(tree);
 xx.Track::Loop(RawRootFile_c);
 cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
 return 1;
 }

