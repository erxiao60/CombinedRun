#define Corr_cxx
#include "Corr.h"
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
#include "math.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"

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
void TrackSeek(Int_t nHits,Int_t *Plane,Int_t *Dimension,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][6],Int_t *PMTdy,Int_t *BGOend, Double_t FitTrack[][9],Double_t FitADC[][8])
{
  int flag=0; 
  int nTlayers=0;
  int nTbars=0;
//  int ngoodFitHits[2]={6,6};
  Float_t MaxADC[2][8];
  Int_t   MaxBar[2][8];//[2]both sides of BGO bar
 for(int o=0;o<2;o++)
for(int p=0;p<12;p++)
{
MaxADC[o][p]=0;
MaxBar[o][p]=0;
} 
 Float_t ADCCutPed=0; 
  for(int i=0; i<nHits;i++)
  for(int j=0;j<2;j++)
  if(BGOend[i]==j)
  {
    nTlayers=Plane[i]*2+Dimension[i];
    nTbars=nTlayers*44+Bar[i]-1;
    flag=(PMTdy[i]-2)/3;
    if(PMTdy[i]==8)//only use Dy8 hits to get a track!!!
    if (ADC[i]>(PedestalPar[nTbars][flag*2]+5*PedestalPar[nTbars][flag*2+1]))
    {
      ADCCutPed=ADC[i]-PedestalPar[nTbars][flag*2];
      if(ADCCutPed>MaxADC[j][nTlayers])  
      {
        MaxADC[j][nTlayers]=ADCCutPed;
        MaxBar[j][nTlayers]=Bar[i]-1;//0--21
      }
// EMC_Mips[nTlayers][Bar[i]][flag]->Fill(ADCCutPed);
    }
  } //single event Loop end
int xy[2]={MaxBar[0][0],MaxBar[0][1]};
FitTrack[0][8]=0;//flag 0:vertical muon 1:squint 
for(int il=0;il<Nlayer;il++)
{
if(MaxBar[0][il]+MaxBar[1][il]==43)
{
if(MaxBar[0][il]==xy[il%2])
{
FitTrack[0][il]=MaxBar[0][il];
FitTrack[1][il]=MaxBar[1][il];
FitADC[0][il]=MaxADC[0][il];
FitADC[1][il]=MaxADC[1][il];
}
//else FitTrack[0][8]=1;
}
else FitTrack[0][8]=1;
}
//if(FitTrack[0][8]==0)
//{
//cout<<"~~~~~~~~~~~~~~~~~~~"<<endl;
//for(int il=0;il<Nlayer;il++)
//for(int is=0;is<2;is++)
//cout<<"Layer:"<<il<<" : "<<MaxBar[is][il]<<"  "<<MaxADC[is][il]<<endl;
//}
//cout<<"vertical muon? "<<FitTrack[0][8];

}


void Corr::Loop(char *RootFileName)
{

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//output .root .eps files
string namesuffix="_Corr.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[50]="./Corr/";
char *EPSFile=strcat(FilePath,outputfile.c_str());
cout<<"Output file: "<<EPSFile<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Define TH1D,TH2D..
TF1 *Corr=new TF1("Corr","landau",0,2500);
TF1 *myCorr=new TF1("myCorr",langaufun,0,2500,4);
Double_t LangauFitpar[Nlayer][Nbar][4];
myCorr->SetParNames("Width","MP","Area","GSigma");

TH1D *CorrCoe[Nlayer][Nbar/2][22];//Correnuation coefficient
for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar/2;ib++)//Bar 0===Bar 44 ,1==43 ,A,B end of a same BGO bar
for(int ip=0;ip<22;ip++)
{
char attc[50];
sprintf(attc,"Layer%d_BGO%d__Point%d",il+1,ib+1,ip);
CorrCoe[il][ib][ip]=new TH1D(attc,attc,300,150,5000);
}

TH2D *Unif[Nlayer][Nbar/2];
for(int il=0;il<Nlayer;il++)
for(int ib=0;ib<Nbar/2;ib++)//Bar 0===Bar 44 ,1==43 ,A,B end of a same BGO bar
{
char coe[50];
sprintf(coe,"Layer%d_BGO%d_Unif",il+1,ib+1);
Unif[il][ib]=new TH2D(coe,coe,25,0,25,220,300,2500);
}

gStyle->SetOptStat(0001);
gStyle->SetOptFit(0111);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Loop start
int BadTrack=0;
int Tevent=0;
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
//void TrackSeek(Int_t nHits,Int_t *Plane,Int_t *Dimension,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][6],Int_t *PMTdy,Int_t *BGOend, Double_t FitTrack[][9],Double_t FitADC[][8])
 Double_t FitTrack[2][9];//FitTrack[][6] :number of well fitted hits
 Double_t FitADC[2][8];
 TrackSeek(nHits,Plane,Dimension,Bar,ADC,PedestalPar,PMTdy,BGOend,FitTrack,FitADC);
//TH1D *CorrCoe[nlayer][Nbar][22];//Correnuation coefficient
int maxbar[2][8];
if(FitTrack[0][8]==0)
for(int ip=0;ip<Nplane;ip++)
{
maxbar[0][ip*2]=(int)FitTrack[0][ip*2];
maxbar[0][ip*2+1]=(int)FitTrack[0][ip*2+1];
CorrCoe[ip*2][maxbar[0][ip*2]][maxbar[0][ip*2+1]]->Fill(sqrt(FitADC[0][ip*2]*FitADC[1][ip*2])); 
CorrCoe[ip*2+1][maxbar[0][ip*2+1]][maxbar[0][ip*2]]->Fill(sqrt(FitADC[0][ip*2+1]*FitADC[1][ip*2+1])); 
}
else BadTrack++;
Tevent++;
}//Loop end

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Corr Fit and Draw
//TH1D *CorrCoe[nlayer][Nbar][22];//Correnuation coefficient
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());

//const int Nlayer=2;//layer = (plane-1)*2+dimension
//const int Nbar=44;//Nbar=0,Nbar=12 for spare channels
//const int Ndy=3;
//const int Nplane=2;
Double_t par[3];
TSpectrum *s=new TSpectrum(1,30);
Int_t nfound =0;
//printf("Found %d candidate peaks \n",nfound);
Float_t *xpeaks;
TCanvas *EMCCorr[Nplane][Ndimension][Nbar/2];
for (int iplane=0;iplane<Nplane;iplane++)
for(int idim=0;idim<Ndimension;idim++)
for(int ibar=0;ibar<Nbar/2;ibar++)
{char att[50];
sprintf(att,"Plane%d_Dim%d_Bar%d",iplane,idim,ibar);
EMCCorr[iplane][idim][ibar]=new TCanvas(att,att,0,0,600,800);
EMCCorr[iplane][idim][ibar]->Divide(4,6);

for(int ip=0;ip<22;ip++)
{
EMCCorr[iplane][idim][ibar]->cd(ip+1);
//~~~~~~~~~~~~~~~~~~~~~~~~~~
nfound=s->Search(CorrCoe[iplane*2+idim][ibar][ip],22,"",0.10);
xpeaks=s->GetPositionX();
if(xpeaks[0]<200);
xpeaks[0]=(float)(CorrCoe[iplane*2+idim][ibar][ip]->GetMean());
Corr->SetRange(xpeaks[0]-300,xpeaks[0]+600);
Corr->SetParameters(100,xpeaks[0]);
CorrCoe[iplane*2+idim][ibar][ip]->Fit(Corr,"R0");
Corr->GetParameters(par);
Corr->SetRange(par[1]-2.5*par[2],par[1]+6*par[2]);
Corr->SetParameters(par);
CorrCoe[iplane*2+idim][ibar][ip]->Fit(Corr,"R0");
Corr->GetParameters(par);
myCorr->SetRange(par[1]-2.0*par[2],par[1]+6.0*par[2]);
//myCorr->SetParameters(par);
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
myCorr->SetParameter(0,par[2]);
myCorr->SetParameter(1,par[1]);
myCorr->SetParameter(3,par[2]);
myCorr->SetLineColor(2);
CorrCoe[iplane*2+idim][ibar][ip]->Draw();
CorrCoe[iplane*2+idim][ibar][ip]->Fit(myCorr,"R");
myCorr->GetParameters(LangauFitpar[iplane*2+idim][ibar]);
//TH2D *Coeff[Nlayer][Nbar];
Unif[iplane*2+idim][ibar]->Fill(ip,LangauFitpar[iplane*2+idim][ibar][1]);
}
EMCCorr[iplane][idim][ibar]->Print(EPSFile);
}

TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Correnuation coefficient
/*TH2D *CoeBar[Nside];//Bar 1=88 (bar+43*layer)
CoeBar[0]=new TH2D("CoeBar0","CoeVsBar_Side0;Tot_Bar;Coe(cm)",176,0,176,50,0,250);
CoeBar[1]=new TH2D("CoeBar1","CoeVsBar_Side1;Tot_Bar;Coe(cm)",176,0,176,50,0,250);
TH1D *CoeDis[Nside];
CoeDis[0]=new TH1D("CoeDis0","CoeDis_Side0;lambda(cm);counts",50,0,250);
CoeDis[1]=new TH1D("CoeDis1","CoeDis_Side1;lambda(cm);counts",50,0,250);
Double_t CorrFunPar[2];
TF1 *CorrFun=new TF1("CorrFun","[0]*exp((37/12)*x/[1])",1,21);//37/12:BGO pitch
*/
TCanvas *ccc[Nlayer];//2 layers, 2 sides
for(int il=0;il<Nlayer;il++)
{char CoeC[30];
sprintf(CoeC,"Layer%d",il+1);
ccc[il]=new TCanvas(CoeC,CoeC,0,0,600,800);
ccc[il]->Divide(4,6);
for(int ib=0;ib<Nbar/2;ib++)
{
ccc[il]->cd(ib+1);
Unif[il][ib]->Draw();
Unif[il][ib]->SetMarkerStyle(23);
Unif[il][ib]->SetMarkerSize(0.3);
Unif[il][ib]->SetMarkerColor(4);
//CorrFun->SetParameters(1200,120*(is*2-1));
//CorrFun->SetLineColor(2);
//Coeff[il][ttbar]->Fit(CorrFun,"R");
//CorrFun->GetParameters(CorrFunPar);
}
if(il==0)
ccc[il]->Print("./Corr/Uniformity.eps(");
if(il==Nlayer-1)
ccc[il]->Print("./Corr/Uniformity.eps)");
else ccc[il]->Print("./Corr/Uniformity.eps");
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~
/*//show coefficient distribution
TCanvas *A=new TCanvas("Coe_Dis","Coe_Dis",0,0,600,600);
A->Divide(2,2);
for(int is=0;is<Nside;is++)
{
A->cd(is*2+1);
CoeBar[is]->Draw();
CoeBar[is]->SetMarkerStyle(23);
CoeBar[is]->SetMarkerSize(0.3);
CoeBar[is]->SetMarkerColor(4);
A->cd(is*2+2);
CoeDis[is]->Draw();
}
A->Print("./Corr/Dis_Coe.eps","Portrait");
*///
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ofstream outputpar;
outputpar.open("./Corr/Corr.txt");
if(!outputpar.good())
{
cout<<"can not open out put fitted par file :./Corr/Corr_LangauPar.txt"<<endl;
exit(1);
}
for(int i=0;i<Nlayer;i++)
for(int j=1;j<Nbar-1;j++)
{
for(int n=0;n<4;n++)
outputpar<<LangauFitpar[i][j][n]<<" ";
outputpar<<"\n";
}
outputpar.close();
cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;

}//Loop end 




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//define main()
string filename;
void helpinfo(){
	cout<<"Usage is ./Corr.exe ./Raw2ROOT/<filename>\n";
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
 
 Corr xx;
 xx.Corr::Init(tree);
 xx.Corr::Loop(RawRootFile_c);
 cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
 return 1;
 }

