#define DyRel_cxx
#include "DyRel.h"
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
void TrackSeek(Int_t nHits,Int_t *Plane,Int_t *Dimension,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][6],Int_t *PMTdy, Double_t FitTrack[][7])
{
  int flag=0; 
  int nTlayers=0;
  int nTbars=0;
  int ngoodFitHits[2]={6,6};
  Float_t MaxADC[12]={'0'};
  Int_t   MaxBar[12]={'0'};
  Float_t ADCCutPed=0; 
  Float_t BeamDir[12]={0,1,2,3,4,5,6,7,8,9,10,11};
  TH2D *Track_z[2];
  Track_z[0]=new TH2D("Track_zx","Track_zx",15,0,15,14,0,14);
  Track_z[1]=new TH2D("Track_zy","Track_zy",15,0,15,14,0,14);
  Double_t LinearPar[2][2];
  for(int i=0; i<nHits;i++)
  {
    nTlayers=Plane[i]*2+Dimension[i];
    nTbars=nTlayers*13+Bar[i];
    flag=(PMTdy[i]-2)/3;
    if(PMTdy[i]==8)//only use Dy8 hits to get a track!!!
    if (ADC[i]>(PedestalPar[nTbars][flag*2]+5*PedestalPar[nTbars][flag*2+1]))
    {
      ADCCutPed=ADC[i]-PedestalPar[nTbars][flag*2];
      if(ADCCutPed>MaxADC[nTlayers])  
      {
        MaxADC[nTlayers]=ADCCutPed;
        MaxBar[nTlayers]=Bar[i];
      }
// EMC_Mips[nTlayers][Bar[i]][flag]->Fill(ADCCutPed);
    }
  } //single event Loop end
  for(int j=0;j<12;j++)
  Track_z[j%2]->Fill(BeamDir[j],MaxBar[j]);
  Track_z[0]->Fit(Tracks);
  Tracks->GetParameters(LinearPar[0]);
  Track_z[1]->Fit(Tracks);
  Tracks->GetParameters(LinearPar[1]);
//  Double_t FitTrack[2][7];//FitTrack[][6] :number of well fitted hits
//  if(LinearPar[0][1]-(int)LinearPar[0][1]-0.5>0.1||LinearPar[0][1]-(int)LinearPar[0][1]-0.5<-0.1)//cut event hit the gaps
//  if(LinearPar[1][1]-(int)LinearPar[1][1]-0.5>0.1||LinearPar[1][1]-(int)LinearPar[1][1]-0.5<-0.1)
//  {
  for(int j=0;j<12;j++)
  {FitTrack[j%2][(j-j%2)/2]=LinearPar[j%2][0]*BeamDir[j]+LinearPar[j%2][1];
  if(FitTrack[j%2][(j-j%2)/2]-MaxBar[j]>1||FitTrack[j%2][(j-j%2)/2]-MaxBar[j]<-1)
  ngoodFitHits[j%2]--;
  }
  FitTrack[0][6]=ngoodFitHits[0];
  FitTrack[1][6]=ngoodFitHits[1];
//  }
delete Track_z[0];
delete Track_z[1];
}



void DyRel::Loop(char *RootFileName)
{

//###########output .root .eps files
string namesuffix="_DyRel.eps";
string outputfile=(string)RootFileName;
outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
int Lname=sizeof(outputfile.c_str());
//Lname-sizeof("./Raw2ROOT/")=Lname-11
string out;
out.assign(outputfile,11,Lname-11);
outputfile=out;
char FilePath[50]="./DyRel/";
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//std::vector<std::vector<TH2D*> > CosDyRel;

  TFile *outFile=new TFile("DyRel.root","RECREATE");
  if (outFile==(TFile*) NULL)
  {
    std::cerr<<"Error:Could not Open Conversion ROOT File!"<<endl;
    exit(1);
  }
  TTree *CosTree= new TTree("AttDy","AttDyCos");  

std::vector<TH2D*>  CosDyRel[Nlayer];
char cc[Nlayer][Nbar][30];
for(Int_t ilayer=0;ilayer<Nlayer;ilayer++)
 { for(Int_t ibar=0;ibar<Nbar;ibar++)
  {
    sprintf(cc[ilayer][ibar],"Layer%d_Bar%d_Dy8/Dy5",(unsigned)(ilayer+1),(unsigned)(ibar+1));
CosDyRel[ilayer].push_back(new TH2D(cc[ilayer][ibar],cc[ilayer][ibar],210,-50,1200,700,-500,16000));
  } 
//CosDyRel.push_back(CosDyRel[ilayer]);
}
gStyle->SetOptStat(111);
gStyle->SetOptFit(0111);
//############# Loop start
int Tevent=0;
Double_t FitRange[Nlayer][Nbar][4];
// Initialize FitRange
for(int ii=0;ii<Nlayer;ii++)
for(int jj=0;jj<Nbar;jj++)
for(int kk=0;kk<4;kk++)
FitRange[ii][jj][kk]=0;
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
//void TrackSeek(Int_t nHits,Int_t *Plane,Int_t *Dimension,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][6],Int_t *PMTdy, Double_t FitTrack[2][7])
/*ouble_t FitTrack[2][7];//FitTrack[][6] :number of well fitted hits
 TrackSeek(nHits,Plane,Dimension,Bar,ADC,PedestalPar,PMTdy,FitTrack);
 
  cout<<"XXXXXXXXXXXXXXXXxx"<<endl;
  cout<<FitTrack[0][6]<<"\n"<<FitTrack[1][6]<<endl;
  cout<<"XXXXXXXXXXXXXXXXxx"<<endl;
  if(FitTrack[0][6]>=4&&FitTrack[1][6]>=4)
   {
    for(int i=0;i<nHits;i++)
      {
      int nTlayers=Plane[i]*2+Dimension[i];
      int nTbars=nTlayers*13+Bar[i];
      int flag=(PMTdy[i]-2)/3;
      if(FitTrack[Dimension[i]][Plane[i]]-Bar[i]>=-0.5&&FitTrack[Dimension[i]][Plane[i]]-Bar[i]<=0.5)
      if (ADC[i]>(PedestalPar[nTbars][flag*2]+5*PedestalPar[nTbars][flag*2+1]))
        {
        Double_t ADCCutPed=ADC[i]-PedestalPar[nTbars][flag*2];
        EMC_Mips[nTlayers][Bar[i]][flag]->Fill(ADCCutPed);
        }
      }
    } //single event Loop end
  else BadTrack++;
  cout<<"XXXXXXXXXXXXXXXXxx"<<endl;
  cout <<BadTrack<<endl;
  cout<<Tevent<<endl;
  cout<<"XXXXXXXXXXXXXXXXxx"<<endl;
  Tevent++;
*/
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
   //Fit Range 
    for(Int_t il=0;il<Nlayer;il++)
    for(Int_t ib=0;ib<Nbar;ib++)
    {
      
    if(Dynode[il][ib][2]>4000&&Dynode[il][ib][2]<5500&&FitRange[il][ib][2]<50)    
    FitRange[il][ib][2]=Dynode[il][ib][1];
    else if(Dynode[il][ib][2]>11000&&Dynode[il][ib][2]<12000&&FitRange[il][ib][3]<50)    
    FitRange[il][ib][3]=Dynode[il][ib][1];
   if(Dynode[il][ib][1]!=0&&Dynode[il][ib][2]!=0)
   {
//if(Dynode[il][ib][2]/Dynode[il][ib][1]>20)
   //cout<<"Dynode:"<<Dynode[il][ib][1]<<"  "<<Dynode[il][ib][2]<<endl;
   CosDyRel[il][ib]->Fill(Dynode[il][ib][1],Dynode[il][ib][2]);
   }
 }   
  for(int il=0;il<Nlayer;il++)
  for(int ib=0;ib<Nbar;ib++)  
  for(int ir=0;ir<Ndy;ir++)
     Dynode[il][ib][ir]=0;
  /*  for(int i=0;i<nHits;i++)
      {
      int nTlayers=Plane[i]*2+Dimension[i];
//      int nTbars=nTlayers*44+Bar[i];
      int flag=(PMTdy[i]-2)/3;
      EMC_Mips[nTlayers][Bar[i]-1][flag]->Fill(ADC[i]);
      }
  */
}//Loop end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//############# DyRel Fit and Draw
TF1 *DyRel=new TF1("DyRel","[0]*x+[1]",50,400);
Double_t Fitpar[Nlayer][Nbar][4];
DyRel->SetParNames("Slope","Intercept");
gStyle->SetOptStat(111);
gStyle->SetOptFit(0111);
TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
string Star=(string)EPSFile;
string FroBracket="[";
Star=Star+FroBracket;
Start->Print(Star.c_str());
//Double_t par[2];
TCanvas *EMCDyRel[Nlayer][Nside];
for (int il=0;il<Nlayer;il++)
for(int is=0;is<Nside;is++)
{
char d[50];
sprintf(d,"Layer%d_Side%d_Dy8/Dy5",il,is);
EMCDyRel[il][is]=new TCanvas(d,d,0,0,600,800);
EMCDyRel[il][is]->Divide(4,6);

for(int ib=0;ib<Nbar/2;ib++)
{
EMCDyRel[il][is]->cd(ib+1);

DyRel->SetRange(FitRange[il][ib+is*22][2],FitRange[il][ib+is*22][3]);
CosDyRel[il][ib+is*22]->Draw();
DyRel->SetLineColor(2);
if (FitRange[il][ib+is*22][2]*FitRange[il][ib+is*22][3]!=0)
  {
  CosDyRel[il][ib+is*22]->Fit(DyRel,"R");
  DyRel->GetParameters(&Fitpar[il][ib+is*22][2]);
  }
  else 
  memset(&Fitpar[il][ib+22*is][2],0,sizeof(Fitpar[il][ib+is*22])/2);
}
EMCDyRel[il][is]->Print(EPSFile);
}


TCanvas *End=new TCanvas("End","End",0,0,600,800);
string En=(string)EPSFile;
string BacBracket="]";
En=En+BacBracket;
End->Print(En.c_str());


ofstream outputpar;
outputpar.open("./DyRel/DyRel_LinearPar.txt");
if(!outputpar.good())
{
cout<<"can not open out put fitted par file :./DyRel/DyRel_LinearPar.txt"<<endl;
exit(1);
}
for(int i=0;i<Nlayer;i++)
for(int j=0;j<Nbar/2;j++)
{
for(int n=2;n<4;n++)
outputpar<<Fitpar[i][j][n]<<" ";
for(int n=2;n<4;n++)
outputpar<<Fitpar[i][43-j][n]<<" ";
outputpar<<"\n";
}
outputpar.close();

cout<<"Total event: "<<Tevent<<endl;

  outFile->cd();
  outFile->Write();
  outFile->Close();
//Loop end 
}



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
 
 DyRel xx;
 xx.DyRel::Init(tree);
 xx.DyRel::Loop(RawRootFile_c);
 cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
 return 1;
 }

