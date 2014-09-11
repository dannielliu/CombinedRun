#define MIPs_cxx
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
TF1 *Tracks=new TF1("Tracks","[0]*x+[1]",0,15); 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//BGO constants information
//Plane(1-7), Layer(1-14), Bar(1-44) ,Side(0,1),Dimension(0,1),FEEcard Number(1-16),FEEchannel Number(144)
//const int Nlayer=7;
//const int Nch=24;
//const int Ndy=3;
//const int Nside=1;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TrackSeek(Int_t nHits,Int_t *Layer,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][2],Int_t *Dy,Int_t *Side,Int_t *GID,Int_t FitTrack[16])
{
  Float_t MaxADC[14];
  Int_t   MaxBar[14];//[2]both sides of BGO bar
  Int_t nH=0;
  memset(MaxADC,0,sizeof(MaxADC));
  memset(MaxBar,0,sizeof(MaxBar));
  Float_t ADCCutPed=0; 
  for(int i=0; i<nHits;i++)//only check side 0
    {
      if(Dy[i]==8&&Side[i]==0)//only use Dy8 hits to get a track!!!
	{
	  if (ADC[i]>(PedestalPar[GID[i]][0]+3*PedestalPar[GID[i]][1]))
	    {
                  nH++;
	      ADCCutPed=ADC[i]-PedestalPar[GID[i]][0];
	      if(ADCCutPed>MaxADC[Layer[i]])  
		{
		  MaxADC[Layer[i]]=ADCCutPed;
		  MaxBar[Layer[i]]=Bar[i];//0--21
		}
	    }
	}
    } //single event Loop end

  if((MaxADC[0]>200||MaxADC[1]>200)||(MaxADC[4]>200||MaxADC[5]>200))
    {
      FitTrack[14]=1;FitTrack[15]=nH;
      for(int il=0;il<Nlayer;il++)
	{
	  FitTrack[il]=MaxBar[il];
	}}
}

void RootConvert::Loop(char *RootFileName){
  gROOT->ProcessLine(".L /home/zhzhy/DRAW/bes3plotstyle.C");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //output .root .eps files
  string namesuffix="_MIPs";
  string outputfile=(string)RootFileName;
  outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
  int Lname=sizeof(outputfile.c_str());
  //Lname-sizeof("./Raw2ROOT/")=Lname-11
  string out;
  out.assign(outputfile,11,Lname-11);
  string TXTFile0="./MIPs/"+out+".txt";
  string EPSFile="./MIPs/"+out+".eps";
  string TXTFile1="./MIPs/MIPs.txt";
  cout<<"Output file 1 "<<EPSFile<<endl;
  cout<<"Output file 2 "<<TXTFile0<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //read pedestal data
  cout<<"Read Pedestal..."<<endl;

  const int nPar=2;//Dy2,5,8 mean and sigma (gaus fitted)
  //const int nGID=(Nlayer*Nside)*Nch*Ndy;
  Double_t PedestalPar[nGID][nPar];
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
  TH1F *EMC_Mips[nGID];
  for (int ng=0;ng<nGID;ng++)
    {
      // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nch+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
      Int_t idy=ng%Ndy;
      Int_t ibar=((int)(ng/Ndy))%Nch;
      Int_t iside=((int)(ng/Nch/Ndy))%Nside;
      Int_t ilayer=(int)(ng/Nside/Nch/Ndy);
      char cc[30];
      sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_MIPs",ilayer,iside,ibar,(idy*3+2));
      if(iside==0||(iside==1&&(ilayer==0||ilayer==13)))
      EMC_Mips[ng]=new TH1F(cc,cc,300,10,8000);
      else
      EMC_Mips[ng]=new TH1F(cc,cc,300,10,4000);
    }

  TH1F *Hits=new TH1F("Hits","Hits;N;counts",616,0,616);
  
 
  TF1 *MIPs=new TF1("MIPs","landau",0,8000);
  TF1 *myMIPs=new TF1("myMIPs",langaufun,0,8000,4);
  myMIPs->SetParNames("Width","MP","Area","GSigma");
  const int nLangaus=4;
  Double_t MIPsPar[nGID][nLangaus];
  for(int ig=0;ig<nGID;ig++)
    {
      memset(MIPsPar[ig],0,sizeof(MIPsPar[ig]));
    }
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(0111);
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
      
      Int_t FitTrack[16];//
      //void TrackSeek(Int_t nHits,Int_t *Layer,Int_t *Bar,Int_t *ADC,Double_t PedestalPar[][2],Int_t *Dy,Int_t *Side,Int_t nGID,Double_t FitTrack[15],Double_t FitADC[14])
      TrackSeek(nHits,Layer,Bar,ADC,PedestalPar,Dy,Side,GID,FitTrack);
      Hits->Fill(FitTrack[15]);
      if(FitTrack[14]==1)
	for(int i=0;i<nHits;i++)
	  {
	    if(FitTrack[Layer[i]]==Bar[i]&&Dy[i]==8)//both sides
	      if (ADC[i]>(PedestalPar[GID[i]][0]+3*PedestalPar[GID[i]][1]))
		EMC_Mips[GID[i]]->Fill(ADC[i]-PedestalPar[GID[i]][0]);
	  }
      else 
	{BadTrack++;
	}
      Tevent++;
    }
  //Loop end
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // MIPs Fit and Draw
  TCanvas *Start=new TCanvas("Start","Start",0,0,600,800);
  string Star=(string)EPSFile;
  string FroBracket="[";
  Star=Star+FroBracket;
  Start->Print(Star.c_str());

  //const int nLangaus=4;
  //Double_t MIPsPar[nGID][nLangaus];
  Double_t par[3];
  //TSpectrum *s=new TSpectrum(1,30);
  //Int_t nfound =0;
  //Float_t *xpeaks;
  TCanvas *EMCMIPs[Nlayer][Nside];
  for(int il=0;il<Nlayer;il++)
    {
      char d[50];
      for(int is=0;is<Nside;is++)
	{
	  sprintf(d,"Layer%d_Side%d",il+1,is);
	  EMCMIPs[il][is]=new TCanvas(d,d,0,0,600,800);
	  EMCMIPs[il][is]->Divide(4,6);
	  for(int ib=0;ib<Nch;ib++)
	    {
	      Int_t iGID=((il*Nside+is)*Nch+ib)*Ndy+2;//idy=2;
	      //mean=EMC_Mips[iplane*2+idim][ibar][0]->GetMean();
	      //RMS =EMC_Mips[iplane*2+idim][ibar][0]->GetRMS();
	      EMCMIPs[il][is]->cd(ib+1);
	      //nfound=s->Search(EMC_Mips[iGID],22,"",0.10);
	      //xpeaks=s->GetPositionX();
	      //if(xpeaks[0]<500);
	      Double_t mean=EMC_Mips[iGID]->GetMean();
              Double_t RMS=EMC_Mips[iGID]->GetRMS();
	      MIPs->SetRange(mean*0.4,mean*1.4);
              //MIPs->SetRange(xpeaks[0]-150,xpeaks[0]+300);
	      MIPs->SetParameters(300,mean*0.8);
	      //MIPs->SetRange(2014-300,2014+600);
	      //MIPs->SetParameters(100,2014);
	      EMC_Mips[iGID]->Fit(MIPs,"R0");
	     // MIPs->GetParameters(par);
	     // MIPs->SetRange(par[1]-3.0*par[2],par[1]+6*par[2]);
	     // MIPs->SetParameters(par);
	     // EMC_Mips[iGID]->Fit(MIPs,"R0");
	      MIPs->GetParameters(par);
	      myMIPs->SetRange(par[1]-3.6*par[2],par[1]+6.50*par[2]);
	      //myMIPs->SetParameters(par);
	      // Once again, here are the Landau * Gaussian parameters:
	      //   par[0]=Width (scale) parameter of Landau density
	      //   par[1]=Most Probable (MP, location) parameter of Landau density
	      //   par[2]=Total area (integral -inf to inf, normalization constant)
	      //   par[3]=Width (sigma) of convoluted Gaussian function
	      myMIPs->SetParameter(0,par[2]);
	      myMIPs->SetParameter(1,par[1]);
	      myMIPs->SetParameter(3,par[2]);
              myMIPs->SetParLimits(3,0,mean);
              myMIPs->SetParLimits(1,mean*0.6,mean);
	      myMIPs->SetLineColor(2);
	      myMIPs->SetLineWidth(2);
	      EMC_Mips[iGID]->Draw();
	      EMC_Mips[iGID]->Fit(myMIPs,"R");
	      myMIPs->GetParameters(MIPsPar[iGID]);
              if(MIPsPar[iGID][1]>mean||MIPsPar[iGID][1]<0.4*mean||MIPsPar[iGID][3]<0||MIPsPar[iGID][3]>MIPsPar[iGID][1]){
	        myMIPs->SetRange(mean*0.4,mean*1.4);
	        myMIPs->SetParameter(0,80);
	        myMIPs->SetParameter(1,mean*0.8);
	        myMIPs->SetParameter(3,300);
	        EMC_Mips[iGID]->Fit(myMIPs,"R");
	        myMIPs->GetParameters(MIPsPar[iGID]);
              }
	    }
	  EMCMIPs[il][is]->Print(EPSFile.c_str());
	}

    }
  TCanvas *End=new TCanvas("End","End",0,0,600,800);
  string En=(string)EPSFile;
  string BacBracket="]";
  En=En+BacBracket;
  End->Print(En.c_str());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ofstream outputpar;
  outputpar.open(TXTFile0.c_str());
  if(!outputpar.good())
    {
      cout<<"can not open out put fitted par file :./MIPs/MIPs_LangauPar.txt"<<endl;
      exit(1);
    }
  for(int i=0;i<nGID;i++)
    {
      if(i%3==2)
	{
	  //const int nLangaus=4;
	  //Double_t MIPsPar[nGID][nLangaus];
	  for(int j=0;j<nLangaus;j++)
	    {
	      outputpar<<MIPsPar[i][j]<<" ";
	    }
	  outputpar<<"\n";
	}
    }
  outputpar.close();
  outputpar.open(TXTFile1.c_str());
  if(!outputpar.good())
    {
      cout<<"can not open out put fitted par file :./MIPs/MIPs_LangauPar.txt"<<endl;
      exit(1);
    }
  for(int i=0;i<nGID;i++)
    {
      if(i%3==2)
	{
	  //const int nLangaus=4;
	  //Double_t MIPsPar[nGID][nLangaus];
	  for(int j=0;j<nLangaus;j++)
	    {
	      outputpar<<MIPsPar[i][j]<<" ";
	    }
	  outputpar<<"\n";
	}
    }
  outputpar.close();
  cout<<"Total event: "<<Tevent<<"\nBad Track event: "<<BadTrack<<endl;
  //Draw nHits
  TCanvas *HH=new TCanvas("HH","HH");
  gPad->SetLogy();
  Hits->Draw();
  HH->Print("nHits/nHits.eps");
  
}//Loop end 




//define main()
string filename;
void helpinfo(){
  cout<<"Usage is ./MIPs.exe ./Raw2ROOT/<filename>\n";
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

