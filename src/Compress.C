#define BeamTest_cxx

#include "vector"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TTree.h"
#include "TH2.h"
#include "TProfile.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include <TFile.h>
#include <string.h>
#include "fstream"
#include "iostream"
#include "Constant.h"
#include <iomanip>
#include <bitset>
using std::bitset;
using std::cout;
using namespace std;
using namespace cosmictest;
string filename;
void helpinfo(){
	cout<<"Usage is ./Compress.exe <filename>\n";
	cout<<"default output file name is <inputfilename.root>"<<endl;	
	return;
}

void phrase_command(int argc,char *argv[]){
	if (argc<2){ helpinfo();
	     exit(0);
	}else {filename=(string) argv[1]; cout<<"Input File is "<<filename<<endl;}
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//main
int main(int argc,char *argv[])
{
  phrase_command(argc,argv);

  // denfine Ntuple
  string namesuffix=".root";
  string outputfile=filename;
  outputfile.replace(strlen(filename.c_str())-4,4,namesuffix);
  cout<<"Output file name is: "<<outputfile<<endl;
  char FilePath[80]="./Raw2ROOT/";
  char *ROOTFile=strcat(FilePath,outputfile.c_str());
  //char *FilePath=strcat("./Raw2ROOT/",outputfile.c_str());
  //TFile *outFile=new TFile(outputfile.c_str(),"RECREATE");
  TFile *outFile=new TFile(ROOTFile,"RECREATE");
  if (outFile==(TFile*) NULL)
  {
    std::cerr<<"Error:Could not Open Conversion ROOT File!"<<endl;
    exit(1);
  }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//BGO constants information
//Nlayer,Nbar,Ndy,Nside,N_FEE;
ifstream BGOinfor;
BGOinfor.open("./map/BGO.infor");
if(!BGOinfor.good())
{
cout<<"~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Can not open BGO constant information file!"<<endl;
exit(-1);
}
int ConBGO[5];
char infN[2][80];
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"~~~~BGOECAL Cosmic Test~~~~"<<endl;
BGOinfor.getline(infN[0],80);
BGOinfor.getline(infN[1],80);
cout<<infN[0]<<endl;
cout<<infN[1]<<endl;
for(int inf=0;inf<5;inf++)
{
BGOinfor>>ConBGO[inf];
cout<<ConBGO[inf]<<"  ";
}
cout<<"\n";
//const int Nlayer=ConBGO[0];
//const int Nbar=ConBGO[1];
//const int Ndy=ConBGO[2];
//const int Nside=ConBGO[3];
//const int N_FEE=ConBGO[4];
//read FEE ID
char infFee[2][80];
Int_t FEEOrder[N_FEE];
BGOinfor.getline(infFee[0],80);
BGOinfor.getline(infFee[1],20);
cout<<"~~~~~~~~FEE cards :~~~~~~~~"<<endl;
for(int fo=0;fo<N_FEE;fo++)
{BGOinfor>>FEEOrder[fo];
cout<<FEEOrder[fo]<<"  ";
}
cout<<"\n";
//read N_channel
char infChan[2][80];
Int_t N_channel[N_FEE];
BGOinfor.getline(infChan[0],80);
BGOinfor.getline(infChan[1],80);
cout<<"~~~~~~~FEE Channel Numbers :~~~~~~~~"<<endl;
for(int fo=0;fo<N_FEE;fo++)
{BGOinfor>>N_channel[fo];
cout<<N_channel[fo]<<"  ";
}
BGOinfor.close();
cout<<"\n";
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//map file
int totalN=504*4;
int n_start[16]={0,144,288,360,504,648,792,864,1008,1152,1296,1368,1512,1656,1800,1872};
//const int UseChan=Nlayer*Nbar*Nside*Ndy; 
  const int UseChan=totalN;
  const int TotalChan=totalN;
  int FEE2EMC[TotalChan][7];
 //initialization
 for(int m=0;m<TotalChan;m++)
 for(int n=0;n<7;n++)
 FEE2EMC[m][n]=0;
  int TChanN=0;
  char mapinfor[80];
  ifstream MapL0;//map of EMC plane&layer ,Dimension,BGO,Side,Dy,FEE_ID ,Chan...
  MapL0.open("./map/map_cosmic");
  if(!MapL0.good())
  {
  cout<<"Open MapL0 file Failed!"<<endl;
  exit(-1);
  }
  MapL0.getline(mapinfor,80);
//  cout<<mapinfor<<endl;
  cout<<"read map infor..."<<endl;
  for(int i=0;i<UseChan;i++)
   {
  MapL0>>TChanN;
  FEE2EMC[TChanN][0]=TChanN;
  for(int j=1;j<7;j++)
  MapL0>>FEE2EMC[TChanN][j];//
   }
  MapL0.close();
  ///////////
  TTree *EMCTree= new TTree("CTest","CTest");  
  Long64_t time=0;
  int event=0;
  int nHits=0;
  int FEE_ID[5000];
  int Chan[5000];
  int ADC[5000];
  int Layer[5000];
  int Bar[5000];//1-22
  int Side[5000];// 1:front or left 2:rear or right
  int Dy[5000];
  int GID[5000];//GID =((layer*2+side)*24+bar)*3;  
  EMCTree->Branch("time",&time,"time/L");
  EMCTree->Branch("event",&event,"event/I");
  EMCTree->Branch("nHits",&nHits,"nHits/I");
  EMCTree->Branch("FEE_ID",FEE_ID,"FEE_ID[nHits]/I");
  EMCTree->Branch("Chan",Chan,"Chan[nHits]/I");
  EMCTree->Branch("ADC",ADC,"ADC[nHits]/I");
  EMCTree->Branch("Layer",Layer,"Layer[nHits]/I");
  EMCTree->Branch("Bar",Bar,"Bar[nHits]/I");
  EMCTree->Branch("Side",Side,"Side[nHits]/I");
  EMCTree->Branch("Dy",Dy,"Dy[nHits]/I");
  EMCTree->Branch("GID",GID,"GID[nHits]/I");

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //define TH1D,TH1F 
  //read data from RawData File
  cout<<"Loop Raw Data..."<<endl;
 //E225 header
  unsigned char header_gl[2];
  unsigned char header_dir[6];
  unsigned char header_time[8];


  Int_t header[4]={0,0,0,0};
  unsigned char  time_b[4]={0};
  Int_t  tail[4]={0,0,0,0};
  Int_t  header_subDAQ[2]={0,0};
  Int_t  tail_subDAQ;
   //for raw data
  signed char data_h;
  unsigned char data_u;
  unsigned char data_l;
   
  Int_t data_2[144];
  memset(data_2,0,sizeof(data_2));
  Int_t ADC_buffer[N_FEE][144]; 
  Int_t Tri_buffer[N_FEE];// FEE triggers buffer
  Int_t chl_buffer[N_FEE][144];
  Int_t hitN[N_FEE];  
  for(int nf=0;nf<N_FEE;nf++) 
  {
  Tri_buffer[nf]=0;
  hitN[nf]=0;
  for(int ch=0;ch<N_channel[nf];ch++) {
  ADC_buffer[nf][ch]=0;
  chl_buffer[nf][ch]=0;
  }
  }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//read binary data
  Int_t Tevent=0;
  Int_t trigger=0;
  Int_t mode=10;//0:normal 1:compress 2:calibration
  bool header_flag=false;
  bool trifill=true;
  bool tailerror=false;
  Int_t SizePac=0;//Size of FEE data packet
  ifstream stream(filename.c_str(), ios_base::binary);
  //ifstream *header_hold;
  for(; !stream.eof();){
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//data check
/*
if(Tevent==54139)
{
for(int ck=0;ck<1000;ck++)
{
for(int ar=0;ar<8;ar++)
{
unsigned char Check;
    for(int k=0;k<2;k++)
    {
    stream.read((char *)(&Check), 1);
Int_t xxx=Check*1;
    cout<<setfill('0')<<setw(2)<<hex<<xxx;
    }
  cout<<" ";
}
cout<<endl;
}
break;
}*/ if(header_flag==false){
    stream.read((char *)(&header_gl[0]),1);
    if(header_gl[0]==0xe2){
    stream.read((char *)(&header_gl[1]),1);
    if(header_gl[1]==0x25){
    for(int dir=0;dir<6;dir++)
    stream.read((char *)(&header_dir[dir]),1);
    for(int tim=0;tim<8;tim++){
    stream.read((char *)(&header_time[tim]),1);
    if(tim==7)
    time=time+header_time[tim];
    else
    time=(time+header_time[tim])*256;
    }
    header_flag=true;
    //cout<<"~~~Tevent~~~"<<Tevent<<endl;
    }else continue;}else continue;}
    
    else{
    stream.read((char *)(&header_subDAQ[0]),1);
    if(header_subDAQ[0]==0xeb){
    stream.read((char *)(&header_subDAQ[1]),1);
    if(header_subDAQ[1]==0xeb){
    stream.read((char *)(&header_subDAQ[1]),1);}
    if(header_subDAQ[1]==0x90){
    
    stream.read((char *)(&header[0]), 1);
    stream.read((char *)(&header[1]), 1);
    stream.read((char *)(&header[2]), 1);
    stream.read((char *)(&header[3]), 1);
//mode :first 2 bits ; FEE ID the other 6 bits
    mode=(Int_t)((Int_t)header[1]/64);//mode=0x01:compressed data; mode=0x00: primary data;
    header[1]=header[1]%64;
//time package check    
    
    SizePac = header[2]*256+header[3]*1;
    //SizePac = SizePac&0x3fff;
    if(SizePac!=294&&SizePac!=150&&mode==0)
    cout<<"Tevent:"<<Tevent<<" FEE"<<header[1]<<" Size:"<<SizePac<<" mode:"<<mode<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//FEE buffer Fill 
    for (int nf=0;nf<N_FEE;nf++)
    {
  //cout<<"~~~~~~~~~~~~~~~~~~~~~"<<endl;
  //  cout<<"FEE:"<<header[1]<<endl;   
      if(header[1]==FEEOrder[nf]){
        if(mode==1){
         int tchl=(int)((SizePac-6)/3);
         if(tchl<=N_channel[nf])//compressed
        for(int iChan=0;iChan<tchl;iChan++) {
          stream.read((char *)(&chl_buffer[nf][hitN[nf]]),1);
          
          stream.read((char *)(&data_h), 1);
         // stream.read((char *)(&data_u), 1);
          stream.read((char *)(&data_l), 1);
          data_2[iChan] = data_h*256+data_l*1;
         // data_2[iChan] = data_u*256+data_l*1;
         // if(data_2[iChan]>32767)
         // data_2[iChan]-=65536;
	  if(chl_buffer[nf][hitN[nf]]<0){
          cout<<"Tevent:"<<Tevent<<" FEE"<<(int)header[1]<<" Size:"<<SizePac<<" mode:"<<mode<<endl;
          cout<<"SizePac:"<<SizePac<<endl;
          cout<<"hitN:"<<hitN[nf]<<endl;      
	  cout<<"chl_buffer:"<<chl_buffer[nf][hitN[nf]]<<endl;}
	  else if(chl_buffer[nf][hitN[nf]]<N_channel[nf]&&chl_buffer[nf][hitN[nf]]>0){
          ADC_buffer[nf][chl_buffer[nf][hitN[nf]]]=data_2[iChan];
          //cout<<data_2[iChan]<<" ";
          hitN[nf]++;
          }
          }
          if((SizePac-6)%3!=0){
          unsigned char spare;
          stream.read((char *)(&spare), 1);
          }
          }
          else if(mode==0||mode==2){//primary
          for(int iChan=0;iChan<N_channel[nf];iChan++) {
          chl_buffer[nf][iChan]=iChan;
          stream.read((char *)(&data_h), 1);
          //stream.read((char *)(&data_u), 1);
          stream.read((char *)(&data_l), 1);
          data_2[iChan] = data_h*256+data_l*1;
          //data_2[iChan] = data_u*256+data_l*1;
          //if(data_2[iChan]>32767)
          //data_2[iChan]-=65536;
          ADC_buffer[nf][iChan]=data_2[iChan];
          hitN[nf]++;
	  }}
          //cout<<"FEE:"<<header[1]<<"  nHit:"<<hitN[header[1]%16]<<endl;
       for(int tai=0;tai<4;tai++) {
         stream.read((char *)(&tail[tai]), 1);
         } 
       Int_t hh=(tail[0]*256+tail[1])%4096;
       hh=hh+1;
    //   if(hh==2960)
     //  cout<<Tevent<<"  "<<Tri_buffer[nf]<<endl;

      if(Tri_buffer[nf]==0)
       Tri_buffer[nf]=hh;
      /* stream.read((char *)(&tail_subDAQ),2);
       if(tail_subDAQ!=42330){
         tailerror=true;
         cout<<"tail error !"<<endl;
      }*/
	  }}
    

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //data loop
    }else continue;//header 0x90 check
     }else continue;//header 0xeb check
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      trifill=true;
      bool trimatch=true;
      for(int tr=0;tr<N_FEE;tr++)
        {
        if(Tri_buffer[tr]==0){
          trifill=false;
          break;
          }
        }
      if(trifill==true)
      for(int tr=0;tr<N_FEE;tr++)
        {
         if(Tri_buffer[tr]!=Tri_buffer[0]){
          cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
          cout<<"Trigger unmatch!"<<endl;
          trimatch=false;
          cout<<"trigger:"<<trigger+1<<" Tevent:"<<Tevent<<endl;
          for(int n=0;n<N_FEE;n++)
            {cout<<Tri_buffer[n]<<"  ";}
          cout<<"\n"<<endl;
          break;
          }
        }
   // trigger cheack 
      
      if(trifill==true&&trimatch==false) {
	//same trigger number
	Int_t tMax=Tri_buffer[0];
	Int_t tMin=Tri_buffer[0];
         for(int n=0;n<N_FEE;n++){
	   if(Tri_buffer[n]>tMax)
             tMax=Tri_buffer[n];
	   if(Tri_buffer[n]<tMin)
             tMin=Tri_buffer[n];
	 }
        if((tMax-tMin)>=2){
        memset(Tri_buffer,0,sizeof(Tri_buffer));
        memset(hitN,0,sizeof(hitN));
	}
        else
        for(int n=0;n<N_FEE;n++){
        if(Tri_buffer[n]!=tMax)
	{Tri_buffer[n]=0;hitN[n]=0;
        }}
      }
      /* int tt=0;
      if(trifill==true&&trimatch==false) {
      trigger=Tri_buffer[0];
      int tri=Tri_buffer[0];
      for(int n=1;n<N_FEE;n++){
      if((Tri_buffer[n]<tri)||(Tri_buffer[n]==4096&&tri==1))
      Tri_buffer[n]=0;
      else if((Tri_buffer[n]>tri&&Tri_buffer[n]<tri+10)||(Tri_buffer[n]==1&&tri==4096)){
      if(Tri_buffer[n]>Tri_buffer[tt])
      trigger=Tri_buffer[n];
      tt=n;}
     else if(Tri_buffer[n]!=tri)
     Tri_buffer[n]=0;
     hitN[n]=0; 
      }
      if(tt!=0){
      memset(Tri_buffer,0,sizeof(Tri_buffer));
      Tri_buffer[tt]=trigger;
         for(int n=0;n<N_FEE;n++)
         {cout<<Tri_buffer[n]<<"  ";}
      }
      }*/    
            
    if(trifill==true&&trimatch==true)
    { 
      if((Tri_buffer[0]-trigger%4096)!=1)
      cout<<"~~~~~"<<"trigger last:"<<trigger<<" trigger now:"<<Tri_buffer[0]<<"~~~~~"<<endl;
      trigger=Tri_buffer[0];
      for(int iFEE=0;iFEE<N_FEE;iFEE++)
      {
      // Int_t nChan=n_start[iFEE];
      Int_t nst=n_start[FEEOrder[iFEE]%16];
      Int_t nChan=0;
      for(int iChan=0;iChan<hitN[iFEE];iChan++)
      {
//	     cout<<hitN[iFEE]<<endl; 
        nChan=nst+chl_buffer[iFEE][iChan];
        if(FEE2EMC[nChan][5]==0)
        continue;
//Tchan   Layer   Side    Bar     Dy      FEE     FeeChannel
        FEE_ID[nHits]=FEE2EMC[nChan][5];
        Chan[nHits]=FEE2EMC[nChan][6];
        ADC[nHits]=ADC_buffer[iFEE][chl_buffer[iFEE][iChan]];
       //cout<<"Event:"<<event<<" "<<hitN[iFEE]<<" "<<ADC[nHits]<<" "<<endl; 
        Layer[nHits]=FEE2EMC[nChan][1];
        Bar[nHits]=FEE2EMC[nChan][3];
        Side[nHits]=FEE2EMC[nChan][2];// 1:front or left 2:rear or right
//      Side[nHits]=0;// 1:front or left 2:rear or right
        Dy[nHits]=FEE2EMC[nChan][4];
        GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
        //GID =((layer*2+side)*24+bar)*3;  
        nHits++; 
        }}
        //cout<<"nHits:"<<nHits<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//reset buffer
       mode=100; 
       for(int iFEE=0;iFEE<N_FEE;iFEE++)
       {
       hitN[iFEE]=0;
       Tri_buffer[iFEE]=0;
       for(int iChan=0;iChan<N_channel[iFEE];iChan++)
       {ADC_buffer[iFEE][iChan]=0;
       chl_buffer[iFEE][iChan]=0;}
       }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Fill Tree  
       EMCTree->Fill();
       event++;
       if(event%1000==0)
       cout<<"eventID is "<<event<<endl;
      //if(event==500) break;
  //Initialization Ntuple
       nHits=0;  
       time=0;	
       memset(FEE_ID,0,sizeof(FEE_ID));
       memset(Chan,0,sizeof(Chan));
       memset(ADC,0,sizeof(ADC));
       memset(Layer,0,sizeof(Layer));
       memset(Bar,0,sizeof(Bar));
       memset(Side,0,sizeof(Side));
       memset(Dy,0,sizeof(Dy));
       memset(GID,0,sizeof(GID));
Tevent=event;
header_flag=false;
  }
    }
}//data read end 
  outFile->cd();
  outFile->Write();
  outFile->Close();
cout<<"Loop end!"<<endl;
cout<<"Total Events:"<<Tevent<<endl;
  return 0;
}
