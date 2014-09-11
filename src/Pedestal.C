#define Pedestal_cxx
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
#include "Pedestal.h"
using namespace std;
using namespace cosmictest;
void RootConvert::Loop(char *RootFileName)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  string namesuffix="_pedestal";
  string outputfile=(string)RootFileName;
  outputfile.replace(strlen(RootFileName)-5,5,namesuffix);
  int Lname=sizeof(outputfile.c_str());
  string out;
  out.assign(outputfile,11,Lname-11);
  string EPSFile="./Pedestal/"+out+".eps";
  string RFile="./Calibration/"+out+".root";
  string RName=out+".root";
  cout<<"Output files:"<<endl;
  cout<<EPSFile<<"\n"<<endl;


  //file name
  ofstream name;
  name.open("./Calibration/filename");
  name<<RName;
  name.close();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //define the TH1F for pedestal 
  TH1F *EMC_Cal[nGID];
  for (int ng=0;ng<nGID;ng++)
    { 
      // GID[nHits]=((Layer[nHits]*Nside+Side[nHits])*Nch+Bar[nHits])*Ndy+(Dy[nHits]-2)/3;
      Int_t idy=ng%Ndy;
      Int_t ibar=((int)(ng/Ndy))%Nch;
      Int_t iside=((int)(ng/Nch/Ndy))%Nside;
      Int_t ilayer=(int)(ng/Nside/Nch/Ndy);
      char cc[30];
      sprintf(cc,"Layer%d_Side%d_Bar%d_Dy%d_Ped",ilayer,iside,ibar,(idy*3+2));
      EMC_Cal[ng]=new TH1F(cc,cc,800,-400,400);
    }
  TF1 *PedFit=new TF1("PedFit","gaus",100,800);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Loop start
  CalPar calpar;
  strcpy(calpar.FileName,(char*)RootFileName);
  CalParIO calpario;
 // calpario.CalParIO::SetPath("./Calibration/Pedestal.root");
  calpario.CalParIO::SetPath(RFile.c_str());
  calpario.CalParIO::InitPedPar(&calpar); 
  
  cout<<"Loop Root File..."<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"~~~~~~~~~~~~~~~"<<nentries<<endl;
  
  Long64_t ientry;
  Long64_t Gstart=0;
  Long64_t Gstop=0;
  Long64_t EntBuf=0;
  Long64_t timegap=600;//10mins
  //cout<<"please input the time duration for pedestal (second):"<<endl;
  //cin>>timegap;
//get Gstart Gstop
  ientry = LoadTree(0);
  fChain->GetEntry(0);
  Gstart=time;

  //for cut last data
  Long64_t tg=Gstart;

 
  ientry = LoadTree(nentries-1);
  fChain->GetEntry(nentries-1);
  Gstop=time;
  cout<<"Gstart & Gstop time: "<<Ast2Dat(Gstart)<<"("<<Gstart<<")  to  "<<Ast2Dat(Gstop)<<"("<<Gstop<<")"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Loop the pedestal par
  Int_t EntNubChk=0;

  Int_t Nloop=(int)((Gstop-Gstart)/timegap)+1;
  
  cout<<"Nloop:"<<Nloop-1<<endl;
  for(int t=0;t<Nloop;t++){ 
    ientry = LoadTree(EntBuf);
    fChain->GetEntry(EntBuf);
    calpar.timestart=time;
     EntNubChk=0;
    //A single loop for pedestal one event
    for (Long64_t jentry=EntBuf; jentry<nentries;jentry++){
      ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      if (mode!=0)  //select data style
      continue;
      
      EntNubChk++;  //check Fill number
 
      //cut last data
      if(time-tg>300){
        Gstart=time;
        cout<<"Cut last package & Reset Gstart:"<<endl;
        cout<<"Gstart:"<<Gstart<<endl;
        calpar.timestart=time;
        for(int iGID=0;iGID<nGID;iGID++){
        EMC_Cal[iGID]->Reset();
      }
      }
      tg=time;


     
      if(time>Gstart+(t+1)*timegap||jentry==nentries-1){
	    EntBuf=jentry;
	    calpar.timestop=time;
            if(jentry==nentries-1){t=Nloop; break; } 
	    if(t<(Nloop-2)&&EntNubChk>=2000)
	      break;
	  }
	
     //loop entries 
      if(jentry%2000==0){
        cout<<"Events:"<<jentry<<endl;
        }
	
     for(int i=0; i<nHits;i++){
       EMC_Cal[GID[i]]->Fill(ADC[i]);
       }
     }
	
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//fitting
cout<<"Fitting..."<<endl;
    Double_t mean=0;
    Double_t par[3];
    //printf("Found %d candidate peaks \n",nfound);
//if(t<Nloop-1)
if(t!=Nloop-3)
for(int iGID=0;iGID<nGID;iGID++){
		  mean=EMC_Cal[iGID]->GetMean();
                  PedFit->SetRange(mean-25,mean+25);
		  PedFit->SetParameters(1000,mean);
		  EMC_Cal[iGID]->Fit(PedFit,"R0Q");
		  PedFit->GetParameters(par);
		  PedFit->SetRange(par[1]-2.5*par[2],par[1]+2.5*par[2]);
		  PedFit->SetParameters(par);
		  PedFit->SetLineColor(2);
		  EMC_Cal[iGID]->Fit(PedFit,"RQ0");
		  //extract the Pedestal parameters
		  PedFit->GetParameters(par);
		  calpar.PedMean[iGID]=par[1];
		  calpar.PedSigma[iGID]=par[2];
                  EMC_Cal[iGID]->Reset();
}
else{
    //Pedestal Fit and Draw
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(111);
    string figure=EPSFile;
    string figure_start=EPSFile+"[";
    string figure_stop=EPSFile+"]";

    TCanvas *figures=new TCanvas("figures","figures",0,0,600,800);
    figures->Print(figure_start.c_str(),"Portrait");

    TCanvas *EMC_Ped[Nlayer][Nside][Ndy];
    for(int ilayer=0;ilayer<Nlayer;ilayer++)
      for(int iside=0;iside<Nside;iside++)
	{
	  char d[3][50];
	  //Dy2,Dy5,Dy8
	  for(int idy=0;idy<Ndy;idy++)
	    {
	      sprintf(d[idy],"Layer%d_Side%d_Dy%d",ilayer+1,iside,idy*3+2);
	      EMC_Ped[ilayer][iside][idy]=new TCanvas(d[idy],d[idy],0,0,600,800);
	      EMC_Ped[ilayer][iside][idy]->Divide(4,6);

	      for(int ibar=0;ibar<Nch;ibar++)
		{
		  Int_t iGID=((ilayer*Nside+iside)*Nch+ibar)*Ndy+idy;
		  EMC_Ped[ilayer][iside][idy]->cd(ibar+1);
		  mean=EMC_Cal[iGID]->GetMean();
                  PedFit->SetRange(mean-25,mean+25);
		  PedFit->SetParameters(1000,mean);
		  EMC_Cal[iGID]->Fit(PedFit,"R0Q");
		  PedFit->GetParameters(par);
		  PedFit->SetRange(par[1]-2.5*par[2],par[1]+2.5*par[2]);
		  PedFit->SetParameters(par);
		  PedFit->SetLineColor(2);
		  EMC_Cal[iGID]->Draw();
		  EMC_Cal[iGID]->Fit(PedFit,"RQ");

		  //extract the Pedestal parameters
		  PedFit->GetParameters(par);
		  calpar.PedMean[iGID]=par[1];
		  calpar.PedSigma[iGID]=par[2];
		  //
 		}
	      EMC_Ped[ilayer][iside][idy]->Print(figure.c_str());
 	    }

 	}
    figures->Print(figure_stop.c_str(),"Portrait");
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   }
    //Set Pedesrtal Par
    //to txt
    char txtfile[60];
    sprintf(txtfile,"./Calibration/Pedestal/P%d-%d",(int)calpar.timestart,(int)calpar.timestop);
    ofstream outtxt;
    outtxt.open(txtfile);
    if(!outtxt.good()){cout<<"Can not open txt file!"<<endl;}
    for(int j=0;j<nGID;j++)
    outtxt<<calpar.PedMean[j]<<" "<<calpar.PedSigma[j]<<endl;
    outtxt.close();

    //to root 
    cout<<"the step time : "<<Ast2Dat(calpar.timestart)<<"("<<calpar.timestart<<")  to  "<<Ast2Dat(calpar.timestop)<<"("<<calpar.timestop<<")"<<endl;
    cout<<"the total time: "<<Ast2Dat(Gstart)<<"("<<Gstart<<")  to  "<<Ast2Dat(Gstop)<<"("<<Gstop<<")"<<endl;
    cout<<"The Event Number: "<<EntNubChk<<"; Rate: "<<(Float_t)EntNubChk/(Float_t)(calpar.timestop-calpar.timestart)<<"Hz"<<endl;  
    calpario.CalParIO::SetPedPar();
  }//Loop end 
    calpario.CalParIO::FinPedPar();
    cout<<"\n~~~~~*****~~~~~\nTHE END"<<endl;
}



//define main()
string filename;
void helpinfo(){
  cout<<"Usage is ./Pedestal.exe ./Raw2ROOT/<filename>\n";
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
  return 1;

} 

