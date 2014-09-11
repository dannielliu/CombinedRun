#define  MipShow_cxx
#include "Mip.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace cosmictest;
void MipShow(){
  //define
  //TFile *MipShow=new TFile("MipShow.root","RECREATE");
  //TTree *MipShwTree=new TTree("MipShwTree","MipShwTree");

  TH2F *ShowMPV[2];
  for(int iside=0;iside<2;iside++){
    char cc[40];
    sprintf(cc,"MPVDrift(%):(Max-Min)/Mid(side%d)",iside);
    char dd[60];
    sprintf(dd,"MPVDrift(%):(Max-Min)/Mid(side%d);Bar;Layer",iside);
    ShowMPV[iside]=new TH2F(cc,dd,22,0,22,14,0,14);
  }
  TH2F *ShowGsigma[2];
  for(int iside=0;iside<2;iside++){
    char cc[40];
    sprintf(cc,"GsigmaDrift(%):(Max-Min)/Mid(side%d)",iside);
    char dd[60];
    sprintf(dd,"GsigmaDrift(%):(Max-Min)/Mid(side%d);Bar;Layer",iside);
    ShowGsigma[iside]=new TH2F(cc,dd,22,0,22,14,0,14);
  }
  TH2F *MipCalMPV[nGbar];
  TH2F *MipCalGsg[nGbar];
  for (int ng=0;ng<nGbar;ng++)
  {
         // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nbar+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
  Int_t ibar=ng%Nbar;
  Int_t iside=((int)(ng/Nbar))%Nside;
  Int_t ilayer=(int)(ng/Nside/Nbar);
  char cc[30];
  char dd[30];
  //for quadrant 1 ,3
  sprintf(cc,"Layer%d_Side%d_Bar%d_MipMPV",ilayer,iside,ibar);
  sprintf(dd,"Layer%d_Side%d_Bar%d_MipGsigma",ilayer,iside,ibar);
  MipCalMPV[ng]=new TH2F(cc,cc,1000,0,12,800,0,8000);
  MipCalMPV[ng]->SetMarkerStyle(24);
  MipCalMPV[ng]->SetMarkerSize(0.13);
  MipCalGsg[ng]=new TH2F(dd,dd,1000,0,12,200,0,2000);
  MipCalGsg[ng]->SetMarkerStyle(25);
  MipCalGsg[ng]->SetMarkerSize(0.13);
  }
  Double_t MPVMax[nGbar];
  Double_t MPVMin[nGbar];
  Double_t GsigmaMax[nGbar];
  Double_t GsigmaMin[nGbar];
  for(int i=0;i<nGbar;i++){
  MPVMax[i]=-10000;MPVMin[i]=10000;
  GsigmaMax[i]=-10000;GsigmaMin[i]=10000;
  }
  //get filename
  ifstream name;
  name.open("filenameM");
  char filename[80];
  name.getline(filename,80);
  name.close();
  //get Mips
  CalPar calmip;
  CalParIO calpar;
  calpar.CalParIO::SetPath(filename);
  calpar.CalParIO::InitMipPar(&calmip);
  Int_t nentries=calpar.fChain->GetEntries();
  Long64_t tb=0;
 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    calpar.fChain->GetEntry(jentry); 
    if(jentry==0)
    tb=calmip.timestart; 
    Float_t timeM=(Float_t)((calmip.timestart+calmip.timestop)/2);
    timeM=(Float_t)(timeM-tb)/3600.;
    for (int ng=0;ng<nGbar;ng++)
    if(calmip.MipMPV[ng]>10000||calmip.MipMPV[ng]<0||calmip.MipGsigma[ng]>2500||calmip.MipGsigma[ng]<0)
      continue;
    else
    {
    MipCalMPV[ng]->Fill(timeM,calmip.MipMPV[ng]);
    MipCalGsg[ng]->Fill(timeM,calmip.MipGsigma[ng]);
    if(calmip.MipMPV[ng]>MPVMax[ng])
    MPVMax[ng]=calmip.MipMPV[ng];
    if(calmip.MipMPV[ng]<MPVMin[ng])
    MPVMin[ng]=calmip.MipMPV[ng];
    if(calmip.MipGsigma[ng]>GsigmaMax[ng])
    GsigmaMax[ng]=calmip.MipGsigma[ng];
    if(calmip.MipGsigma[ng]<GsigmaMin[ng])
    GsigmaMin[ng]=calmip.MipGsigma[ng];
//  MipCalMPV[ng]->Fill(Ast2Dat((calmip.timestop+calmip.timestart)/2),calmip.MipMPV[ng]);
//  MipCalGsg[ng]->Fill(Ast2Dat((calmip.timestop+calmip.timestart)/2),calmip.MipGsigma[ng]);
    }
  }
  calpar.dataIO->Close();
  
  Int_t BadBuf[nGbar];
  Int_t nBad=0;
  for (int ng=0;ng<nGbar;ng++){
    int il=(int)(ng/Nside/Nbar);
    int ib=(int)ng%Nbar;
    int is=(int)(ng/Nbar)%Nside;
   //cout<<MPVMax[ng]<<"  "<<MPVMin[ng]<<endl;
   // cout<<GsigmaMax[ng]<<"  "<<GsigmaMin[ng]<<endl;
    MipCalGsg[ng]->SetAxisRange((GsigmaMax[ng]+GsigmaMin[ng])/2-500,(GsigmaMax[ng]+GsigmaMin[ng])/2+500,"Y");
    MipCalMPV[ng]->SetAxisRange((MPVMax[ng]+MPVMin[ng])/2-500,(MPVMax[ng]+MPVMin[ng])/2+500,"Y");
    Double_t MPVdrift=(MPVMax[ng]-MPVMin[ng])/((MPVMax[ng]+MPVMin[ng])/2)*100;
    Double_t Gsigmadrift=(GsigmaMax[ng]-GsigmaMin[ng])/((GsigmaMax[ng]+GsigmaMin[ng])/2)*100;    
    ShowMPV[is]->SetBinContent(ib+1,il+1,(Float_t)((int)(MPVdrift*10))/10);
    ShowGsigma[is]->SetBinContent(ib+1,il+1,(Float_t)((int)(Gsigmadrift*10))/10);
    //if(MPVdrift>=10||Gsigmadrift>=10){
    if(MPVdrift>=5||Gsigmadrift>=50){
    BadBuf[nBad]=ng;
    nBad++;
    }
  }
  gStyle->SetOptStat(000);
 TCanvas *Mip_inf0=new TCanvas("Mips inf0","Mips inf0",0,0,800,400);
//  gStyle->SetTextFont(72);
  gStyle->SetTextSizePixels(12);
  Mip_inf0->Divide(2,1);
  for(int iside=0;iside<2;iside++){
      Mip_inf0->cd(iside+1);
      gPad->SetGrid();
      ShowMPV[iside]->Draw("ColTEXTZ");
  }
 Mip_inf0->Print("MipShow.eps(","Portrait");
 TCanvas *Mip_inf1=new TCanvas("Mips inf1","Mips inf1",0,0,800,400);
//  gStyle->SetTextFont(72);
  gStyle->SetTextSizePixels(12);
  Mip_inf1->Divide(2,1);
  for(int iside=0;iside<2;iside++){
      Mip_inf1->cd(iside+1);
      gPad->SetGrid();
      ShowGsigma[iside]->Draw("ColTEXTZ");
  }

 //badShow
 if(nBad==0)
   Mip_inf1->Print("MipShow.eps)");
 else{
   Mip_inf1->Print("MipShow.eps");
   TCanvas *showMS=new TCanvas("showMS","showMS",0,0,800,600);
   showMS->Divide(4,3);
   cout<<"unstable channel number:"<<nBad<<endl;  
   if(nBad>6){ 
     cout<<"fisrt 6 shown in MipShow.eps!!!"<<endl;
     nBad=6;
   }

   char xtitle[50];
   char *begin=Ast2Dat(tb);
   sprintf(xtitle,"time (hours from %s)",begin);

 
   for(int i=0;i<nBad;i++){
     cout<<BadBuf[i]<<endl;
     showMS->cd(2*i+1);
     MipCalMPV[BadBuf[i]]->Draw();
     MipCalMPV[BadBuf[i]]->SetXTitle(xtitle);
     MipCalMPV[BadBuf[i]]->SetYTitle("Mips MPV (ADC channels)");
     
     showMS->cd(2*i+2);
     MipCalGsg[BadBuf[i]]->Draw();
     MipCalGsg[BadBuf[i]]->SetXTitle(xtitle);
     MipCalGsg[BadBuf[i]]->SetYTitle("Mips Gsigma (ADC channels)");
   }
 showMS->Print("MipShow.eps)");
 }




 //MipShwTree->Fill();
 // MipShow->cd();
 // MipShow->Write();
 // MipShow->Close();
}
