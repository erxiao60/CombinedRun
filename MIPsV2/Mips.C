
#include "../include/Constant.h"
using namespace std;
using namespace cosmictest;
int Mips(){
  //define
  TH1D *MIPs_dis[nGside];
  TH2D *MIPs_uni[nGside];
  for(int il=0;il<Nlayer;il++){
    for(int is=0;is<Nside;is++){
      char m[30];
      sprintf(m,"Layer%d_Side%d",il+1,is);
      MIPs_dis[il*Nside+is]=new TH1D(m,m,55,0,1100);
      sprintf(m,"Layer%d_Side%d_unif",il+1,is);
      MIPs_uni[il*Nside+is]=new TH2D(m,m,Nbar,0,Nbar,160,0,2);
    }}
  //read MIPs and Fill
  const int nMPar=4;//Dy2,5,8 mean and sigma (gaus fitted)
  Double_t MIPsPar[nGbar][nMPar];
  ifstream MipsPar;
  MipsPar.open("MIPsV2.txt");
  if (!MipsPar.good())
    {cout<<"Can not open MIPs TXT File!!!"<<endl;
      exit(0);
    }
  //Double_t MIPsPar[nTBar][nDyPar];
  Double_t Max[nGside];
  memset(Max,0,sizeof(Max));
  Double_t Min[nGside];
  for(int ii=0;ii<nGside;ii++)
    Min[ii]=10000.;
  for(int ig=0;ig<nGbar;ig++){
    MipsPar>>MIPsPar[ig][0]>>MIPsPar[ig][1]>>MIPsPar[ig][2]>>MIPsPar[ig][3];
    //cout<<MIPsPar[ig][0]<<MIPsPar[ig][1]<<MIPsPar[ig][2]<<MIPsPar[ig][3]<<endl;
    int ib=ig%Nbar;
    int iGside=(int)(ig/Nbar);
    if(ib!=0&&ib!=23){
      MIPs_dis[iGside]->Fill(MIPsPar[ig][1]);
      if(MIPsPar[ig][1]>Max[iGside])
	Max[iGside]=MIPsPar[ig][1];
      if(MIPsPar[ig][1]<Min[iGside])
	Min[iGside]=MIPsPar[ig][1];
    }
  }
  MipsPar.close();
  for(int ig=0;ig<nGbar;ig++){
    int ib=ig%Nbar;
    int iGside=(int)(ig/Nbar);
    if(ib!=0&&ib!=23){
      //cout<<"Max:"<<Max[iGside]<<" "<<"Min:"<<Min[iGside]<<endl;
      Double_t unif=MIPsPar[ig][1]/(Max[iGside]+Min[iGside])*2;
      MIPs_uni[iGside]->Fill(ib,unif);
    }
  }
  //Draw

  TCanvas *side[Nside][2];//for odd-layer and even-layer 
  for(int i=0;i<Nside;i++){
    for(int j=0;j<2;j++){
      char ss[30];
      //char latex[30];
      if(j==0){
	if(i==0)
	  sprintf(ss,"Quadrant-%d",4);
	if(i==1)
	  sprintf(ss,"Quadrant-%d",2);
      }
      else{
	if(i==0)
	  sprintf(ss,"Quadrant-%d",1);
	if(i==1)
	  sprintf(ss,"Quadrant-%d",3);
      }
      side[i][j]=new TCanvas(ss,ss);
      side[i][j]->Divide(4,4);
      for(int k=0;k<Nlayer;k++){
	if(k%2!=j){
	  side[i][j]->cd(k+j);
	  MIPs_dis[k*Nside+i]->Draw();
	  MIPs_dis[k*Nside+i]->SetLineColor(kBlue);
	  side[i][j]->cd(k+j+1);
	  MIPs_uni[k*Nside+i]->SetMarkerStyle(24);
	  MIPs_uni[k*Nside+i]->SetMarkerSize(0.3);
	  MIPs_uni[k*Nside+i]->SetMarkerColor(kMagenta);
	  MIPs_uni[k*Nside+i]->Draw();
	}
      }
      side[i][j]->cd(15);
      TLegend *l=new TLegend(0.1,0.3,0.9,0.7);
      l->SetTextFont(72);
      l->SetHeader(ss);
      l->SetTextSize(0.20);
      l->SetTextColor(4);
      l->SetFillColor(kYellow-9);
      l->Draw();
    }}
  return 222;
}
