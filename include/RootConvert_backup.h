//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  8 11:19:01 2013 by ROOT version 5.28/00c
// from TTree VA32Data/EMC_VA32Data
// found on file: test_15-27-44_2013-04-12 1000V_01.root
//////////////////////////////////////////////////////////
#ifndef RootConvert_h
#define RootConvert_h

#include "Langaus.h"
#include "Constant.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include <CalParMem.h>
class RootConvert {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Long64_t           time;
   Short_t         mtime;
   Int_t           nHits;
   Int_t          trigger;
   Int_t            mode;
   Int_t           FEE_ID[5000];   //[nHits]
   Int_t           Chan[5000];   //[nHits]
   Int_t           ADC[5000];   //[nHits]
   Int_t           Layer[5000];   //[nHits]
   Int_t           GID[5000];
   Int_t           Bar[5000];
   Int_t           Side[5000];
   Int_t           Dy[5000];

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_time;   //!
   TBranch        *b_mtime;
   TBranch        *b_trigger;
   TBranch        *b_nHits;   //!
   TBranch        *b_mode;
   TBranch        *b_FEE_ID;   //!
   TBranch        *b_Chan;   //!
   TBranch        *b_ADC;   //!
   TBranch        *b_Layer;   //!
   TBranch        *b_GID;   //!
   TBranch        *b_Bar;   //!
   TBranch        *b_Side;   //!
   TBranch        *b_Dy;   //!

   RootConvert();
   virtual ~RootConvert();
  // virtual void     RootPath(TTree *tree=0,Text_t *RootFileName);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(char *RootFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     RawTrack(CalPar &calpar,Int_t *Lmax);
};


RootConvert::RootConvert()
{}
RootConvert::~RootConvert()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootConvert::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootConvert::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RootConvert::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("mtime", &mtime, &b_mtime);
   fChain->SetBranchAddress("mode", &mode, &b_mode);  
   fChain->SetBranchAddress("nHits", &nHits, &b_nHits);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("FEE_ID", FEE_ID, &b_FEE_ID);
   fChain->SetBranchAddress("Chan", Chan, &b_Chan);
   fChain->SetBranchAddress("ADC", ADC, &b_ADC);
   fChain->SetBranchAddress("Layer", Layer, &b_Layer);
   fChain->SetBranchAddress("GID", GID, &b_GID);
   fChain->SetBranchAddress("Bar", Bar, &b_Bar);
   fChain->SetBranchAddress("Side", Side, &b_Side);
   fChain->SetBranchAddress("Dy", Dy, &b_Dy);
   Notify();
}

Bool_t RootConvert::Notify()
{
 return kTRUE;
}

void RootConvert::Show(Long64_t entry)
{
if (!fChain) return;
fChain->Show(entry);
}
Int_t RootConvert::Cut(Long64_t entry)
{
return 1;
}

void RootConvert::RawTrack(CalPar &calpar,Int_t *Lmax){
//  Lmax=(Int_t*)malloc(20);
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
	  if (ADC[i]>calpar.PedMean[GID[i]]+3*calpar.PedSigma[GID[i]]) {
              nH++;
	      ADCCutPed=ADC[i]-calpar.PedMean[GID[i]];
	      if(ADCCutPed>MaxADC[Layer[i]])  
		{
		  MaxADC[Layer[i]]=ADCCutPed;
		  MaxBar[Layer[i]]=Bar[i];//0--21
		}
	    }
	}
    } //single event Loop end
  if((MaxADC[0]>200||MaxADC[1]>200)&&(MaxADC[12]>200||MaxADC[13]>200)) {
  Lmax[14]=1;Lmax[15]=nH;
    for(int il=0;il<Nlayer;il++){
      Lmax[il]=MaxBar[il];
    }
  }
}
#endif // #ifdef RootConvert_cxx
