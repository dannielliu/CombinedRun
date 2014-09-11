#include <TString.h>
#include "Mip.h"
using namespace std;
using namespace cosmictest;
int Mip(){
  //define
  TH1F *MPV_dis=new TH1F("Mips_MPV","MPV Distribution;mean;counts",550,0,5500);
  TH2F *MPV_chan=new TH2F("MipsMPV_vs_Channels","MPV vs Chan;channels(Gbar);mean",nGbar,0,nGbar,550,0,5500);
  TH2F *MPV_layer=new TH2F("","",nGbar,0,nGbar,550,0,5500);
  for (int i=1;i<Nlayer;i++){
    for(int j=0;j<550;j++){
      MPV_layer->Fill(i*Nside*Nbar,j*10);
    }}

  TH1F *Gsigma_dis=new TH1F("Mips_sigma","Gsigma Distribution;sigma;counts",250,0,2500);
  TH2F *Gsigma_chan=new TH2F("MipsGsigma_vs_Channels","Gsigma vs Chan;channels(Gbar);sigma",nGbar,0,nGbar,250,0,2500);
  TH2F *Gsigma_layer=new TH2F("","",nGbar,0,nGbar,250,0,2500);
  for (int i=1;i<Nlayer;i++){
    for(int j=0;j<250;j++){
      Gsigma_layer->Fill(i*Nside*Nbar,j*10);
    }}
  TH2F *Show[2];
  for(int iside=0;iside<2;iside++){
    char cc[40];
    sprintf(cc,"MipsGsigma(side%d_dy%d)",iside,8);
    char dd[60];
    sprintf(dd,"MipsGsigma(side%d_dy%d);Bar;Layer",iside,8);
    Show[iside]=new TH2F(cc,dd,22,0,22,14,0,14);
  }

  //read pedestal and Fill
  const int nDyPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
  //const int nGbar=(Nlayer*Nside)*Nbar*Ndy;
  Double_t MipsPar[nGbar][nDyPar];
//  ifstream MipPar;
//:  MipPar.open("Mips.txt");
//  if (!MipPar.good())
//    {cout<<"Can not open Mips TXT File!!!"<<endl;
//      exit(0);
//    }
  //get filename
  ifstream name;
  name.open("filenameM");
  char file[50];
  name.getline(file,50);
  name.close();
CalPar calmip;
CalParIO calpar;
calpar.CalParIO::SetPath(file);
calpar.CalParIO::InitMipPar(&calmip);
char date[40]="2014-06-17-20-00-00";
//cout<<"Please input the time of the pedestal(YYYY-MM-DD-HH-MM-SS):"<<endl;
//cin>>date; 
calmip.runtime=Dat2Ast(date);

calpar.CalParIO::GetMipPar(&calmip);
calpar.dataIO->Close();

char *begin=Ast2Dat(calmip.timestart);
char *end=Ast2Dat(calmip.timestop);
cout<<"The start and end time of this Mips:"<<endl;
cout<<"Begin:"<<begin<<"("<<calmip.timestart<<")"<<endl;
cout<<"End  :"<<end<<"("<<calmip.timestop<<")"<<endl;
char filename[100];
sprintf(filename,"%s_to_%s",begin,end);

  //Double_t MipsPar[nTBar][nDyPar];
  for(int ig=0;ig<nGbar;ig++){
//  MipPar>>MipsPar[ig][0]>>MipsPar[ig][1];
MipsPar[ig][0]=calmip.MipMPV[ig];
MipsPar[ig][1]=calmip.MipGsigma[ig];
    int il=(int)(ig/Nside/Nbar);
    int ib=ig%Nbar;
    int is=(int)(ig/Nbar)%Nside;
      MPV_dis->Fill(MipsPar[ig][0]);
      MPV_chan->Fill(ig,MipsPar[ig][0]);
      Gsigma_dis->Fill(MipsPar[ig][1]);
      Gsigma_chan->Fill(ig,MipsPar[ig][1]);
/*    if(MipsPar[ig][1]>=20){
      cout<<"~~~~~Mips noise~~~~~~"<<endl;
      cout<<"Mips Gsigma:"<<MipsPar[ig][1]<<endl;
      cout<<"Layer :"<<(int)(ig/Nside/Nbar/Ndy)<<endl;
      cout<<"Side  :"<<(int)(ig/Nbar/Ndy)%Nside<<endl;
      cout<<"Bar   :"<<(int)(ig/Ndy)%Nbar<<endl;
      cout<<"Dy    :"<<ig%Ndy<<endl;
    }
*/
    Show[is]->SetBinContent(ib+1,il+1,(Float_t)((int)(MipsPar[ig][0]*1))/1);
  }

//  MipPar.close();
  //read file name

  //Draw()
  gStyle->SetOptStat(000);
  
 TCanvas *Mip_inf0=new TCanvas("Mips inf0","Mips inf0",0,0,800,400);
//  gStyle->SetTextFont(72);
  gStyle->SetTextSizePixels(12);
  Mip_inf0->Divide(2,1);
  for(int iside=0;iside<2;iside++){
      Mip_inf0->cd(iside+1);
      gPad->SetGrid();
      Show[iside]->Draw("ColTEXTZ");
  }
  TCanvas *Mip_inf=new TCanvas("Mips inf","Mips inf");
  Mip_inf->Divide(2,2);
  Mip_inf->cd(1);
  MPV_dis->Draw();
  MPV_dis->SetLineColor(kMagenta+1);
  TLegend *leg1=new TLegend(0.1,0.80,0.65,0.90);
  leg1->SetTextFont(72);
  leg1->SetHeader(filename);
  leg1->SetTextSize(0.04);
  leg1->SetFillColor(kYellow-9);
  leg1->Draw();
  //gStyle->SetOptStat(111);
  Mip_inf->cd(2);
  MPV_chan->SetMarkerStyle(25);
  MPV_chan->SetMarkerColor(kMagenta+1);
  MPV_chan->SetMarkerSize(0.3);
  MPV_chan->Draw();
  MPV_layer->SetMarkerStyle(23);
  MPV_layer->SetMarkerColor(kBlue);
  MPV_layer->SetMarkerSize(0.2);
  MPV_layer->Draw("SAME");
  TLegend *leg2=new TLegend(0.6,0.65,0.88,0.85);
  leg2->SetTextFont(72);
  leg2->SetHeader("Mips MPV(gaus)");
  leg2->SetTextSize(0.04);
  leg2->SetFillColor(kYellow-9);
  leg2->AddEntry(Gsigma_chan,"Used channels","p");
  leg2->AddEntry(Gsigma_layer,"Layer fences","p");
  leg2->Draw();

  Mip_inf->cd(3);
  Gsigma_dis->Draw();
  Gsigma_dis->SetLineColor(kMagenta+1);
  Float_t SigM=(Float_t)Gsigma_dis->GetMean();
  char SigT[40];
  sprintf(SigT,"UsedChannel (MPV:%.3f)",SigM);

  TLegend *l=new TLegend(0.5,0.65,0.88,0.85);
  l->SetTextFont(72);
  l->SetHeader("Mips Gsigma(gaus)");
  l->SetTextSize(0.04);
  l->SetFillColor(kYellow-9);
  l->AddEntry(Gsigma_dis,SigT,"lp");
  l->Draw();
  gPad->SetLogy();

  Mip_inf->cd(4);
  Gsigma_chan->SetMarkerStyle(25);
  Gsigma_chan->SetMarkerColor(kMagenta+1);
  Gsigma_chan->SetMarkerSize(0.3);
  Gsigma_chan->Draw();
  Gsigma_layer->SetMarkerStyle(23);
  Gsigma_layer->SetMarkerColor(kBlue);
  Gsigma_layer->SetMarkerSize(0.2);
  Gsigma_layer->Draw("SAME");
  TLegend *leg=new TLegend(0.6,0.65,0.88,0.85);
  leg->SetTextFont(72);
  leg->SetHeader("Mips Gsigma(gaus)");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kYellow-9);
  leg->AddEntry(Gsigma_chan,"Used channels","p");
  leg->AddEntry(Gsigma_layer,"Layer fences","p");
  leg->Draw();

  return 222;
}
