//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 29 18:28:55 2016 by ROOT version 5.34/26
// from TTree h9010/
// found on file: ./hmsoptics/hms68734.1.root
//////////////////////////////////////////////////////////

#ifndef HMS_tree_h
#define HMS_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class HMS_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         hcer_npe;
   Float_t         hsp;
   Float_t         hse;
   Float_t         charge;
   Float_t         hsdelta;
   Float_t         hstheta;
   Float_t         hsphi;
   Float_t         w;
   Float_t         hszbeam;
   Float_t         hsdedx1;
   Float_t         hsbeta;
   Float_t         hsshtrk;
   Float_t         hsprtrk;
   Float_t         hsxfp;
   Float_t         hsyfp;
   Float_t         hsxpfp;
   Float_t         hsypfp;
   Float_t         hsytar;
   Float_t         hsxptar;
   Float_t         hsyptar;
   Float_t         hstart;
   Float_t         hsfptime;
   Float_t         eventid;
   Float_t         ev_type;
   Float_t         trigtype;
   Float_t         S0x1padc;
   Float_t         S0x1nadc;
   Float_t         S0x2padc;
   Float_t         S0x2nadc;
   Float_t         S0x1ptdc;
   Float_t         S0x1ntdc;
   Float_t         S0x2ptdc;
   Float_t         S0x2ntdc;
   Float_t         rast_y;
   Float_t         rast_x;
   Float_t         hdchits1;
   Float_t         hdchits2;
   Float_t         hntrkfp;
   Float_t         ntrkhits;
   Float_t         chi2ndf;
   Float_t         fry_cm;
   Float_t         frx_cm;
   Float_t         xtar;

   // List of branches
   TBranch        *b_hcer_npe;   //!
   TBranch        *b_hsp;   //!
   TBranch        *b_hse;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_hsdelta;   //!
   TBranch        *b_hstheta;   //!
   TBranch        *b_hsphi;   //!
   TBranch        *b_w;   //!
   TBranch        *b_hszbeam;   //!
   TBranch        *b_hsdedx1;   //!
   TBranch        *b_hsbeta;   //!
   TBranch        *b_hsshtrk;   //!
   TBranch        *b_hsprtrk;   //!
   TBranch        *b_hsxfp;   //!
   TBranch        *b_hsyfp;   //!
   TBranch        *b_hsxpfp;   //!
   TBranch        *b_hsypfp;   //!
   TBranch        *b_hsytar;   //!
   TBranch        *b_hsxptar;   //!
   TBranch        *b_hsyptar;   //!
   TBranch        *b_hstart;   //!
   TBranch        *b_hsfptime;   //!
   TBranch        *b_eventid;   //!
   TBranch        *b_ev_type;   //!
   TBranch        *b_trigtype;   //!
   TBranch        *b_S0x1padc;   //!
   TBranch        *b_S0x1nadc;   //!
   TBranch        *b_S0x2padc;   //!
   TBranch        *b_S0x2nadc;   //!
   TBranch        *b_S0x1ptdc;   //!
   TBranch        *b_S0x1ntdc;   //!
   TBranch        *b_S0x2ptdc;   //!
   TBranch        *b_S0x2ntdc;   //!
   TBranch        *b_rast_y;   //!
   TBranch        *b_rast_x;   //!
   TBranch        *b_hdchits1;   //!
   TBranch        *b_hdchits2;   //!
   TBranch        *b_hntrkfp;   //!
   TBranch        *b_ntrkhits;   //!
   TBranch        *b_chi2ndf;   //!
   TBranch        *b_fry_cm;   //!
   TBranch        *b_frx_cm;   //!
   TBranch        *b_xtar;   //!

   HMS_tree(TTree *tree=0);
   virtual ~HMS_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HMS_tree_cxx
HMS_tree::HMS_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./hmsoptics/hms68734.1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./hmsoptics/hms68734.1.root");
      }
      f->GetObject("h9010",tree);

   }
   Init(tree);
}

HMS_tree::~HMS_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HMS_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HMS_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HMS_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("hcer_npe", &hcer_npe, &b_hcer_npe);
   fChain->SetBranchAddress("hsp", &hsp, &b_hsp);
   fChain->SetBranchAddress("hse", &hse, &b_hse);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("hsdelta", &hsdelta, &b_hsdelta);
   fChain->SetBranchAddress("hstheta", &hstheta, &b_hstheta);
   fChain->SetBranchAddress("hsphi", &hsphi, &b_hsphi);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("hszbeam", &hszbeam, &b_hszbeam);
   fChain->SetBranchAddress("hsdedx1", &hsdedx1, &b_hsdedx1);
   fChain->SetBranchAddress("hsbeta", &hsbeta, &b_hsbeta);
   fChain->SetBranchAddress("hsshtrk", &hsshtrk, &b_hsshtrk);
   fChain->SetBranchAddress("hsprtrk", &hsprtrk, &b_hsprtrk);
   fChain->SetBranchAddress("hsxfp", &hsxfp, &b_hsxfp);
   fChain->SetBranchAddress("hsyfp", &hsyfp, &b_hsyfp);
   fChain->SetBranchAddress("hsxpfp", &hsxpfp, &b_hsxpfp);
   fChain->SetBranchAddress("hsypfp", &hsypfp, &b_hsypfp);
   fChain->SetBranchAddress("hsytar", &hsytar, &b_hsytar);
   fChain->SetBranchAddress("hsxptar", &hsxptar, &b_hsxptar);
   fChain->SetBranchAddress("hsyptar", &hsyptar, &b_hsyptar);
   fChain->SetBranchAddress("hstart", &hstart, &b_hstart);
   fChain->SetBranchAddress("hsfptime", &hsfptime, &b_hsfptime);
   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("ev_type", &ev_type, &b_ev_type);
   fChain->SetBranchAddress("trigtype", &trigtype, &b_trigtype);
   fChain->SetBranchAddress("S0x1padc", &S0x1padc, &b_S0x1padc);
   fChain->SetBranchAddress("S0x1nadc", &S0x1nadc, &b_S0x1nadc);
   fChain->SetBranchAddress("S0x2padc", &S0x2padc, &b_S0x2padc);
   fChain->SetBranchAddress("S0x2nadc", &S0x2nadc, &b_S0x2nadc);
   fChain->SetBranchAddress("S0x1ptdc", &S0x1ptdc, &b_S0x1ptdc);
   fChain->SetBranchAddress("S0x1ntdc", &S0x1ntdc, &b_S0x1ntdc);
   fChain->SetBranchAddress("S0x2ptdc", &S0x2ptdc, &b_S0x2ptdc);
   fChain->SetBranchAddress("S0x2ntdc", &S0x2ntdc, &b_S0x2ntdc);
   fChain->SetBranchAddress("rast_y", &rast_y, &b_rast_y);
   fChain->SetBranchAddress("rast_x", &rast_x, &b_rast_x);
   fChain->SetBranchAddress("hdchits1", &hdchits1, &b_hdchits1);
   fChain->SetBranchAddress("hdchits2", &hdchits2, &b_hdchits2);
   fChain->SetBranchAddress("hntrkfp", &hntrkfp, &b_hntrkfp);
   fChain->SetBranchAddress("ntrkhits", &ntrkhits, &b_ntrkhits);
   fChain->SetBranchAddress("chi2ndf", &chi2ndf, &b_chi2ndf);
   fChain->SetBranchAddress("fry_cm", &fry_cm, &b_fry_cm);
   fChain->SetBranchAddress("frx_cm", &frx_cm, &b_frx_cm);
   fChain->SetBranchAddress("xtar", &xtar, &b_xtar);
   Notify();
}

Bool_t HMS_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HMS_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HMS_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HMS_tree_cxx
