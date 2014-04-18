#define myadd_cxx
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <utility>
#include <stdio.h>
#include <iostream>
#include <TSystem.h>
#include <vector>
#include "TInterpreter.h"
#include <TStyle.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
using namespace std;

char outfile[50];
char file1[50];
char file2[50];
void helpinfo(){
	cout<<"Usage is ./myadd.exe outfile added-file1 added-file2\n";
//	cout<<"default output file name is <inputfilename.eps>"<<endl;	
	return;
}
void phrase_command(int argc,char *argv[]){
	if (argc<2){ helpinfo();
	     exit(0);
	}else {sprintf(outfile,"%s",argv[1]);
              sprintf(file1,"%s",argv[2]);
              sprintf(file2,"%s",argv[3]);
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"output file: "<<outfile<<endl;
cout<<"added file 1: "<<file1<<endl;
cout<<"added file 2: "<<file2<<endl;
}
}
int main(int argc,char *argv[]){
  phrase_command(argc,argv);
//read root

  TFile *outroot=new TFile(outfile,"RECREATE");
  if (outroot==(TFile*) NULL)
  {
    std::cerr<<"Error:Could not Open merging ROOT File!"<<endl;
    exit(1);
  }
  TTree *EMCTree= new TTree("CTest","CTest");  
  int event[2]={0,0};
  int nHits[2]={0,0};
  int FEE_ID[2][5000];
  int Chan[2][5000];
  int ADC[2][5000];
  int Layer[2][5000];
  int Bar[2][5000];//1-22
  int Side[2][5000];// 1:front or left 2:rear or right
  int Dy[2][5000];
  int GID[2][5000];//GID =((layer*2+side)*24+bar)*3;  
  EMCTree->Branch("event",&event[0],"event/I");
  EMCTree->Branch("nHits",&nHits[0],"nHits/I");
  EMCTree->Branch("FEE_ID",FEE_ID[0],"FEE_ID[nHits]/I");
  EMCTree->Branch("Chan",Chan[0],"Chan[nHits]/I");
  EMCTree->Branch("ADC",ADC[0],"ADC[nHits]/I");
  EMCTree->Branch("Layer",Layer[0],"Layer[nHits]/I");
  EMCTree->Branch("Bar",Bar[0],"Bar[nHits]/I");
  EMCTree->Branch("Side",Side[0],"Side[nHits]/I");
  EMCTree->Branch("Dy",Dy[0],"Dy[nHits]/I");
  EMCTree->Branch("GID",GID[0],"GID[nHits]/I");
//set file 01
   TFile *input_f1= 0;
   if (!gSystem->AccessPathName(file1)) {
      input_f1= TFile::Open(file1);
   } else {
     cout<<"could not find the input file1 !! "<<endl;
      exit(-1); 
   }
 if (!input_f1) return 1;
   TTree *fChain_f1 = (TTree *) input_f1->Get("CTest");
   fChain_f1->SetBranchAddress("event", &event[0] );
   fChain_f1->SetBranchAddress("nHits", &nHits[0]);
   fChain_f1->SetBranchAddress("FEE_ID",FEE_ID[0]);
   fChain_f1->SetBranchAddress("Chan", Chan[0]);
   fChain_f1->SetBranchAddress("ADC", ADC[0]);
   fChain_f1->SetBranchAddress("Layer", Layer[0]);
   fChain_f1->SetBranchAddress("GID", GID[0]);
   fChain_f1->SetBranchAddress("Bar", Bar[0]);
   fChain_f1->SetBranchAddress("Side",Side[0]);
   fChain_f1->SetBranchAddress("Dy", Dy[0]);
//set file 02     
   TFile *input_f2= 0;
   if (!gSystem->AccessPathName(file2)) {
      input_f2 = TFile::Open(file2);
   } else {
     cout<<"could not find the input file2 !! "<<endl;
      exit(-1); 
   }
 if (!input_f2) return 1;
   TTree *fChain_f2 = (TTree *) input_f2->Get("CTest");
   fChain_f2->SetBranchAddress("event", &event[1] );
   fChain_f2->SetBranchAddress("nHits", &nHits[1]);
   fChain_f2->SetBranchAddress("FEE_ID",FEE_ID[1]);
   fChain_f2->SetBranchAddress("Chan", Chan[1]);
   fChain_f2->SetBranchAddress("ADC", ADC[1]);
   fChain_f2->SetBranchAddress("Layer", Layer[1]);
   fChain_f2->SetBranchAddress("GID", GID[1]);
   fChain_f2->SetBranchAddress("Bar", Bar[1]);
   fChain_f2->SetBranchAddress("Side",Side[1]);
   fChain_f2->SetBranchAddress("Dy", Dy[1]);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//merge the event
Int_t entries[2];
entries[0]=fChain_f1->GetEntries();
entries[1]=fChain_f2->GetEntries();
if(entries[0]!=entries[1]){
cout<<"Trigger numbers are not matched!"<<endl;
cout<<"file1:"<<entries[0]<<"  "<<"file2:"<<entries[1]<<endl;
exit(-1);
}
for(int i=0;i<entries[0];i++){
fChain_f1->GetEntry();
fChain_f2->GetEntry();
for(int j=nHits[0];j<nHits[0]+nHits[1];j++){
int k=j-nHits[0];
FEE_ID[0][j]=FEE_ID[1][k];
Chan[0][j]=Chan[1][k];
ADC[0][j]=ADC[1][k];
Layer[0][j]=Layer[1][k];
GID[0][j]=GID[1][k];
Bar[0][j]=Bar[1][k];
Dy[0][j]=Dy[1][k];
Side[0][j]=Side[1][k];
}
nHits[0]=nHits[0]+nHits[1];

EMCTree->Fill();
}


outroot->cd();
outroot->Write();
outroot->Close();
  return 1;
}



