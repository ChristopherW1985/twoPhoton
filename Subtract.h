//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 10 10:39:49 2013 by ROOT version 5.34/00
// from TTree tree/tree
// found on file: twoPhoton.root
//////////////////////////////////////////////////////////

#ifndef Subtract_h
#define Subtract_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Subtract {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        e1;
   Double_t        e2;
   Double_t        tdiff;
   Int_t           pairNo;
   Int_t           veto;

   // List of branches
   TBranch        *b_e1;   //!
   TBranch        *b_e2;   //!
   TBranch        *b_tdiff;   //!
   TBranch        *b_pairNo;   //!
   TBranch        *b_veto;   //!

   Subtract(TTree *tree=0);
   virtual ~Subtract();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     FillHisto();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t checkSequVeto(const int angleGr, const double e);
   virtual Bool_t checkCascade(const double e);
  
};

#endif

#ifdef Subtract_cxx
Subtract::Subtract(TTree *tree) : fChain(0) 
{
   if (tree == 0) 
   {
      TChain * chain = new TChain("tree","MEGA Tree");
      chain->Add("twoPhoton_r206.root/tree");
      chain->Add("twoPhoton_r207.root/tree");
      chain->Add("twoPhoton_r208.root/tree");
      chain->Add("twoPhoton_r214.root/tree");
      chain->Add("twoPhoton_r215.root/tree");
      chain->Add("twoPhoton_r216.root/tree");
      chain->Add("twoPhoton_r217.root/tree");
      chain->Add("twoPhoton_r221.root/tree");
      chain->Add("twoPhoton_r222.root/tree");
      chain->Add("twoPhoton_r223.root/tree");
      chain->Add("twoPhoton_r224.root/tree");
      chain->Add("twoPhoton_r242.root/tree");
      chain->Add("twoPhoton_r254.root/tree");

      /*chain->Add("twoPhoton.root/tree");
	chain->Add("twoPhoton_old.root/tree");
	chain->Add("twoPhoton_new.root/tree");*/
      tree = chain;
   }
   Init(tree);
}

Subtract::~Subtract()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Subtract::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Subtract::LoadTree(Long64_t entry)
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

void Subtract::Init(TTree *tree)
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

   fChain->SetBranchAddress("e1", &e1, &b_e1);
   fChain->SetBranchAddress("e2", &e2, &b_e2);
   fChain->SetBranchAddress("tdiff", &tdiff, &b_tdiff);
   fChain->SetBranchAddress("pairNo", &pairNo, &b_pairNo);
   fChain->SetBranchAddress("veto", &veto, &b_veto);
   Notify();
}

Bool_t Subtract::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Subtract::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Subtract::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Subtract_cxx
