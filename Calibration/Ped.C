#include <TString.h>
#include "Ped.h"
using namespace std;
using namespace cosmictest;
int Ped(){
  //define
  TH1F *Mean_dis=new TH1F("Pedestal_Mean","Mean Distribution;mean;counts",400,-400,400);
  TH1F *Mean_spare=new TH1F("","Spare",400,-400,400);
  TH2F *Mean_chan=new TH2F("PedestalMean_vs_Channels","Mean vs Chan;channels(GID);mean",nGID,0,nGID,400,-400,400);
  TH2F *Mean_sc=new  TH2F("","",nGID,0,nGID,400,-400,400);
  TH2F *Mean_layer=new TH2F("","",nGID,0,nGID,400,-400,400);
  for (int i=1;i<Nlayer;i++){
    for(int j=-50;j<50;j++){
      Mean_layer->Fill(i*Nside*Nch*Ndy,j*8);
    }}

  TH1F *Sigma_dis=new TH1F("Pedestal_sigma","Sigma Distribution;sigma;counts",100,0,20);
  TH1F *Sigma_spare=new TH1F("","Spare",100,0,20);
  TH2F *Sigma_chan=new TH2F("PedestalSigma_vs_Channels","Sigma vs Chan;channels(GID);sigma",nGID,0,nGID,100,0,20);
  TH2F *Sigma_sc=new  TH2F("","",nGID,0,nGID,100,0,20);
  TH2F *Sigma_layer=new TH2F("","",nGID,0,nGID,100,0,20);
  for (int i=1;i<Nlayer;i++){
    for(int j=0;j<100;j++){
      Sigma_layer->Fill(i*Nside*Nch*Ndy,j*0.2);
    }}
  TH2F *Show[2][3];
  for(int iside=0;iside<2;iside++){
    for(int idy=0;idy<3;idy++){
    char cc[40];
    sprintf(cc,"PedestalSigma(side%d_dy%d)",iside,idy*3+2);
    char dd[60];
    sprintf(dd,"PedestalSigma(side%d_dy%d);Bar;Layer",iside,idy*3+2);
    Show[iside][idy]=new TH2F(cc,dd,24,0,24,14,0,14);
    }
  }

  //read pedestal and Fill
  const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
  //const int nGID=(Nlayer*Nside)*Nch*Ndy;
  Double_t PedestalPar[nGID][nDyPar];
//  ifstream PedPar;
//:  PedPar.open("Pedestal.txt");
//  if (!PedPar.good())
//    {cout<<"Can not open Pedestal TXT File!!!"<<endl;
//      exit(0);
//    }
  //get filename
  ifstream name;
  name.open("filename");
  char filename[50];
  name.getline(filename,50);
  name.close();
CalPar calped;
CalParIO calpar;
calpar.CalParIO::SetPath(filename);
calpar.CalParIO::InitPedPar(&calped);
char date[40]="2014-06-17-20-00-00";
//cout<<"Please input the time of the pedestal(YYYY-MM-DD-HH-MM-SS):"<<endl;
//cin>>date; 
calped.runtime=Dat2Ast(date);
calpar.CalParIO::GetPedPar(&calped);
char *begin=Ast2Dat(calped.timestart);
char *end=Ast2Dat(calped.timestop);
cout<<"The start and end time of this Pedestal:"<<endl;
cout<<"Begin:"<<begin<<"("<<calped.timestart<<")"<<endl;
cout<<"End  :"<<end<<"("<<calped.timestop<<")"<<endl;
char filename[100];
sprintf(filename,"%s_to_%s",begin,end);

  //Double_t PedestalPar[nTBar][nDyPar];
  for(int ig=0;ig<nGID;ig++){
//  PedPar>>PedestalPar[ig][0]>>PedestalPar[ig][1];
PedestalPar[ig][0]=calped.PedMean[ig];
PedestalPar[ig][1]=calped.PedSigma[ig];
    int il=(int)(ig/Nside/Nch/Ndy);
    int ib=(int)(ig/Ndy)%Nch;
    int is=(int)(ig/Nch/Ndy)%Nside;
    int id=ig%3;
    if(ib==0||ib==23){
      Mean_spare->Fill(PedestalPar[ig][0]);
      Mean_sc->Fill(ig,PedestalPar[ig][0]);
      Sigma_spare->Fill(PedestalPar[ig][1]);
      Sigma_sc->Fill(ig,PedestalPar[ig][1]);
    }
    else{
      Mean_dis->Fill(PedestalPar[ig][0]);
      Mean_chan->Fill(ig,PedestalPar[ig][0]);
      Sigma_dis->Fill(PedestalPar[ig][1]);
      Sigma_chan->Fill(ig,PedestalPar[ig][1]);
    }
    if(PedestalPar[ig][1]>=20){
      cout<<"~~~~~Pedestal noise~~~~~~"<<endl;
      cout<<"Pedestal Sigma:"<<PedestalPar[ig][1]<<endl;
      cout<<"Layer :"<<(int)(ig/Nside/Nch/Ndy)<<endl;
      cout<<"Side  :"<<(int)(ig/Nch/Ndy)%Nside<<endl;
      cout<<"Bar   :"<<(int)(ig/Ndy)%Nch<<endl;
      cout<<"Dy    :"<<ig%Ndy<<endl;
    }
    Show[is][id]->SetBinContent(ib+1,il+1,(Float_t)((int)(PedestalPar[ig][1]*10))/10);
  }

//  PedPar.close();
  //read file name

  //Draw()
  gStyle->SetOptStat(000);
  
 TCanvas *Ped_inf0=new TCanvas("Pedestal inf0","Pedestal inf0",0,0,1200,800);
//  gStyle->SetTextFont(72);
  gStyle->SetTextSizePixels(12);
  Ped_inf0->Divide(3,2);
  for(int iside=0;iside<2;iside++){
    for(int idy=0;idy<3;idy++){
      Ped_inf0->cd(iside*3+idy+1);
      gPad->SetGrid();
      Show[iside][idy]->Draw("ColTEXTZ");
      }
  }
  TCanvas *Ped_inf=new TCanvas("Pedestal inf","Pedestal inf");
  Ped_inf->Divide(2,2);
  Ped_inf->cd(1);
  Mean_dis->Draw();
  Mean_dis->SetLineColor(kMagenta+1);
  Mean_spare->Draw("SAME");
  Mean_spare->SetLineColor(kGreen+4);
  TLegend *leg1=new TLegend(0.1,0.80,0.65,0.90);
  leg1->SetTextFont(72);
  leg1->SetHeader(filename);
  leg1->SetTextSize(0.04);
  leg1->SetFillColor(kYellow-9);
  leg1->Draw();
  //gStyle->SetOptStat(111);
  Ped_inf->cd(2);
  Mean_chan->SetMarkerStyle(25);
  Mean_chan->SetMarkerColor(kMagenta+1);
  Mean_chan->SetMarkerSize(0.3);
  Mean_chan->Draw();
  Mean_sc->SetMarkerStyle(26);
  Mean_sc->SetMarkerColor(kGreen+4);
  Mean_sc->SetMarkerSize(0.3);
  Mean_sc->Draw("SAME");
  Mean_layer->SetMarkerStyle(23);
  Mean_layer->SetMarkerColor(kBlue);
  Mean_layer->SetMarkerSize(0.2);
  Mean_layer->Draw("SAME");
  TLegend *leg2=new TLegend(0.6,0.65,0.88,0.85);
  leg2->SetTextFont(72);
  leg2->SetHeader("Pedestal Mean(gaus)");
  leg2->SetTextSize(0.04);
  leg2->SetFillColor(kYellow-9);
  leg2->AddEntry(Sigma_chan,"Used channels","p");
  leg2->AddEntry(Sigma_sc,"Spare channels","p");
  leg2->AddEntry(Sigma_layer,"Layer fences","p");
  leg2->Draw();

  Ped_inf->cd(3);
  Sigma_dis->Draw();
  Sigma_dis->SetLineColor(kMagenta+1);
  Sigma_spare->Draw("SAME");
  Sigma_spare->SetLineColor(kGreen+4);
  Float_t SigM=(Float_t)Sigma_dis->GetMean();
  char SigT[40];
  sprintf(SigT,"UsedChannel (Mean:%.3f)",SigM);
  Float_t SpaM=(Float_t)Sigma_spare->GetMean();
  char SpaT[40];
  sprintf(SpaT,"UsedChannel (Mean:%.3f)",SpaM);

  TLegend *l=new TLegend(0.5,0.65,0.88,0.85);
  l->SetTextFont(72);
  l->SetHeader("Pedestal Sigma(gaus)");
  l->SetTextSize(0.04);
  l->SetFillColor(kYellow-9);
  l->AddEntry(Sigma_dis,SigT,"lp");
  l->AddEntry(Sigma_spare,SpaT,"lp");
  l->Draw();
  gPad->SetLogy();

  Ped_inf->cd(4);
  Sigma_chan->SetMarkerStyle(25);
  Sigma_chan->SetMarkerColor(kMagenta+1);
  Sigma_chan->SetMarkerSize(2);
  Sigma_chan->Draw();
  Sigma_sc->SetMarkerStyle(26);
  Sigma_sc->SetMarkerColor(kGreen+4);
  Sigma_sc->SetMarkerSize(0.3);
  Sigma_sc->Draw("SAME");
  Sigma_layer->SetMarkerStyle(23);
  Sigma_layer->SetMarkerColor(kBlue);
  Sigma_layer->SetMarkerSize(0.2);
  Sigma_layer->Draw("SAME");
  TLegend *leg=new TLegend(0.6,0.65,0.88,0.85);
  leg->SetTextFont(72);
  leg->SetHeader("Pedestal Sigma(gaus)");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kYellow-9);
  leg->AddEntry(Sigma_chan,"Used channels","p");
  leg->AddEntry(Sigma_sc,"Spare channels","p");
  leg->AddEntry(Sigma_layer,"Layer fences","p");
  leg->Draw();

  return 222;
}
