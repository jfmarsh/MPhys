//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 13 13:36:01 2021 by ROOT version 6.22/02
// from TTree events/Events tree
// found on file: output_allCalo_reco_noise_40GeV_14layers.root
//////////////////////////////////////////////////////////

#ifndef FourteenLayers40GeV_h
#define FourteenLayers40GeV_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class FourteenLayers40GeV {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxmergedECalCells = 2237;
   static constexpr Int_t kMaxemptyCaloCells = 1;
   static constexpr Int_t kMaxmergedECalCellsNoise = 1716352;
   static constexpr Int_t kMaxCaloClusters = 7;
   static constexpr Int_t kMaxCaloClusters0 = 1;
   static constexpr Int_t kMaxECalBarrelCells = 2237;
   static constexpr Int_t kMaxGenParticles = 1;
   static constexpr Int_t kMaxGenParticles0 = 1;
   static constexpr Int_t kMaxGenParticles1 = 1;

   // Declaration of leaf types
   Int_t           mergedECalCells_;
   ULong64_t       mergedECalCells_core_cellId[kMaxmergedECalCells];   //[mergedECalCells_]
   Float_t         mergedECalCells_core_energy[kMaxmergedECalCells];   //[mergedECalCells_]
   Float_t         mergedECalCells_core_time[kMaxmergedECalCells];   //[mergedECalCells_]
   UInt_t          mergedECalCells_core_bits[kMaxmergedECalCells];   //[mergedECalCells_]
   Int_t           emptyCaloCells_;
   ULong64_t       emptyCaloCells_core_cellId[kMaxemptyCaloCells];   //[emptyCaloCells_]
   Float_t         emptyCaloCells_core_energy[kMaxemptyCaloCells];   //[emptyCaloCells_]
   Float_t         emptyCaloCells_core_time[kMaxemptyCaloCells];   //[emptyCaloCells_]
   UInt_t          emptyCaloCells_core_bits[kMaxemptyCaloCells];   //[emptyCaloCells_]
   Int_t           mergedECalCellsNoise_;
   ULong64_t       mergedECalCellsNoise_core_cellId[kMaxmergedECalCellsNoise];   //[mergedECalCellsNoise_]
   Float_t         mergedECalCellsNoise_core_energy[kMaxmergedECalCellsNoise];   //[mergedECalCellsNoise_]
   Float_t         mergedECalCellsNoise_core_time[kMaxmergedECalCellsNoise];   //[mergedECalCellsNoise_]
   UInt_t          mergedECalCellsNoise_core_bits[kMaxmergedECalCellsNoise];   //[mergedECalCellsNoise_]
   Int_t           CaloClusters_;
   Float_t         CaloClusters_core_energy[kMaxCaloClusters];   //[CaloClusters_]
   Float_t         CaloClusters_core_time[kMaxCaloClusters];   //[CaloClusters_]
   Float_t         CaloClusters_core_position_x[kMaxCaloClusters];   //[CaloClusters_]
   Float_t         CaloClusters_core_position_y[kMaxCaloClusters];   //[CaloClusters_]
   Float_t         CaloClusters_core_position_z[kMaxCaloClusters];   //[CaloClusters_]
   UInt_t          CaloClusters_core_bits[kMaxCaloClusters];   //[CaloClusters_]
   UInt_t          CaloClusters_hits_begin[kMaxCaloClusters];   //[CaloClusters_]
   UInt_t          CaloClusters_hits_end[kMaxCaloClusters];   //[CaloClusters_]
   Int_t           CaloClusters0_;
   Int_t           CaloClusters0_index[kMaxCaloClusters0];   //[CaloClusters#0_]
   Int_t           CaloClusters0_collectionID[kMaxCaloClusters0];   //[CaloClusters#0_]
   Int_t           ECalBarrelCells_;
   ULong64_t       ECalBarrelCells_core_cellId[kMaxECalBarrelCells];   //[ECalBarrelCells_]
   Float_t         ECalBarrelCells_core_energy[kMaxECalBarrelCells];   //[ECalBarrelCells_]
   Float_t         ECalBarrelCells_core_time[kMaxECalBarrelCells];   //[ECalBarrelCells_]
   UInt_t          ECalBarrelCells_core_bits[kMaxECalBarrelCells];   //[ECalBarrelCells_]
   Int_t           GenParticles_;
   Int_t           GenParticles_core_pdgId[kMaxGenParticles];   //[GenParticles_]
   Int_t           GenParticles_core_charge[kMaxGenParticles];   //[GenParticles_]
   UInt_t          GenParticles_core_status[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_vertex_x[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_vertex_y[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_vertex_z[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_p4_mass[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_p4_px[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_p4_py[kMaxGenParticles];   //[GenParticles_]
   Float_t         GenParticles_core_p4_pz[kMaxGenParticles];   //[GenParticles_]
   UInt_t          GenParticles_core_bits[kMaxGenParticles];   //[GenParticles_]
   Int_t           GenParticles0_;
   Int_t           GenParticles0_index[kMaxGenParticles0];   //[GenParticles#0_]
   Int_t           GenParticles0_collectionID[kMaxGenParticles0];   //[GenParticles#0_]
   Int_t           GenParticles1_;
   Int_t           GenParticles1_index[kMaxGenParticles1];   //[GenParticles#1_]
   Int_t           GenParticles1_collectionID[kMaxGenParticles1];   //[GenParticles#1_]

   // List of branches
   TBranch        *b_mergedECalCells_;   //!
   TBranch        *b_mergedECalCells_core_cellId;   //!
   TBranch        *b_mergedECalCells_core_energy;   //!
   TBranch        *b_mergedECalCells_core_time;   //!
   TBranch        *b_mergedECalCells_core_bits;   //!
   TBranch        *b_emptyCaloCells_;   //!
   TBranch        *b_emptyCaloCells_core_cellId;   //!
   TBranch        *b_emptyCaloCells_core_energy;   //!
   TBranch        *b_emptyCaloCells_core_time;   //!
   TBranch        *b_emptyCaloCells_core_bits;   //!
   TBranch        *b_mergedECalCellsNoise_;   //!
   TBranch        *b_mergedECalCellsNoise_core_cellId;   //!
   TBranch        *b_mergedECalCellsNoise_core_energy;   //!
   TBranch        *b_mergedECalCellsNoise_core_time;   //!
   TBranch        *b_mergedECalCellsNoise_core_bits;   //!
   TBranch        *b_CaloClusters_;   //!
   TBranch        *b_CaloClusters_core_energy;   //!
   TBranch        *b_CaloClusters_core_time;   //!
   TBranch        *b_CaloClusters_core_position_x;   //!
   TBranch        *b_CaloClusters_core_position_y;   //!
   TBranch        *b_CaloClusters_core_position_z;   //!
   TBranch        *b_CaloClusters_core_bits;   //!
   TBranch        *b_CaloClusters_hits_begin;   //!
   TBranch        *b_CaloClusters_hits_end;   //!
   TBranch        *b_CaloClusters0_;   //!
   TBranch        *b_CaloClusters0_index;   //!
   TBranch        *b_CaloClusters0_collectionID;   //!
   TBranch        *b_ECalBarrelCells_;   //!
   TBranch        *b_ECalBarrelCells_core_cellId;   //!
   TBranch        *b_ECalBarrelCells_core_energy;   //!
   TBranch        *b_ECalBarrelCells_core_time;   //!
   TBranch        *b_ECalBarrelCells_core_bits;   //!
   TBranch        *b_GenParticles_;   //!
   TBranch        *b_GenParticles_core_pdgId;   //!
   TBranch        *b_GenParticles_core_charge;   //!
   TBranch        *b_GenParticles_core_status;   //!
   TBranch        *b_GenParticles_core_vertex_x;   //!
   TBranch        *b_GenParticles_core_vertex_y;   //!
   TBranch        *b_GenParticles_core_vertex_z;   //!
   TBranch        *b_GenParticles_core_p4_mass;   //!
   TBranch        *b_GenParticles_core_p4_px;   //!
   TBranch        *b_GenParticles_core_p4_py;   //!
   TBranch        *b_GenParticles_core_p4_pz;   //!
   TBranch        *b_GenParticles_core_bits;   //!
   TBranch        *b_GenParticles0_;   //!
   TBranch        *b_GenParticles0_index;   //!
   TBranch        *b_GenParticles0_collectionID;   //!
   TBranch        *b_GenParticles1_;   //!
   TBranch        *b_GenParticles1_index;   //!
   TBranch        *b_GenParticles1_collectionID;   //!

   FourteenLayers40GeV(TTree *tree=0);
   virtual ~FourteenLayers40GeV();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FourteenLayers40GeV_cxx
FourteenLayers40GeV::FourteenLayers40GeV(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_allCalo_reco_noise_40GeV_14layers.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_allCalo_reco_noise_40GeV_14layers.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

FourteenLayers40GeV::~FourteenLayers40GeV()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FourteenLayers40GeV::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FourteenLayers40GeV::LoadTree(Long64_t entry)
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

void FourteenLayers40GeV::Init(TTree *tree)
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

   fChain->SetBranchAddress("mergedECalCells", &mergedECalCells_, &b_mergedECalCells_);
   fChain->SetBranchAddress("mergedECalCells.core.cellId", mergedECalCells_core_cellId, &b_mergedECalCells_core_cellId);
   fChain->SetBranchAddress("mergedECalCells.core.energy", mergedECalCells_core_energy, &b_mergedECalCells_core_energy);
   fChain->SetBranchAddress("mergedECalCells.core.time", mergedECalCells_core_time, &b_mergedECalCells_core_time);
   fChain->SetBranchAddress("mergedECalCells.core.bits", mergedECalCells_core_bits, &b_mergedECalCells_core_bits);
   fChain->SetBranchAddress("emptyCaloCells", &emptyCaloCells_, &b_emptyCaloCells_);
   fChain->SetBranchAddress("emptyCaloCells.core.cellId", &emptyCaloCells_core_cellId, &b_emptyCaloCells_core_cellId);
   fChain->SetBranchAddress("emptyCaloCells.core.energy", &emptyCaloCells_core_energy, &b_emptyCaloCells_core_energy);
   fChain->SetBranchAddress("emptyCaloCells.core.time", &emptyCaloCells_core_time, &b_emptyCaloCells_core_time);
   fChain->SetBranchAddress("emptyCaloCells.core.bits", &emptyCaloCells_core_bits, &b_emptyCaloCells_core_bits);
   fChain->SetBranchAddress("mergedECalCellsNoise", &mergedECalCellsNoise_, &b_mergedECalCellsNoise_);
   fChain->SetBranchAddress("mergedECalCellsNoise.core.cellId", mergedECalCellsNoise_core_cellId, &b_mergedECalCellsNoise_core_cellId);
   fChain->SetBranchAddress("mergedECalCellsNoise.core.energy", mergedECalCellsNoise_core_energy, &b_mergedECalCellsNoise_core_energy);
   fChain->SetBranchAddress("mergedECalCellsNoise.core.time", mergedECalCellsNoise_core_time, &b_mergedECalCellsNoise_core_time);
   fChain->SetBranchAddress("mergedECalCellsNoise.core.bits", mergedECalCellsNoise_core_bits, &b_mergedECalCellsNoise_core_bits);
   fChain->SetBranchAddress("CaloClusters", &CaloClusters_, &b_CaloClusters_);
   fChain->SetBranchAddress("CaloClusters.core.energy", CaloClusters_core_energy, &b_CaloClusters_core_energy);
   fChain->SetBranchAddress("CaloClusters.core.time", CaloClusters_core_time, &b_CaloClusters_core_time);
   fChain->SetBranchAddress("CaloClusters.core.position.x", CaloClusters_core_position_x, &b_CaloClusters_core_position_x);
   fChain->SetBranchAddress("CaloClusters.core.position.y", CaloClusters_core_position_y, &b_CaloClusters_core_position_y);
   fChain->SetBranchAddress("CaloClusters.core.position.z", CaloClusters_core_position_z, &b_CaloClusters_core_position_z);
   fChain->SetBranchAddress("CaloClusters.core.bits", CaloClusters_core_bits, &b_CaloClusters_core_bits);
   fChain->SetBranchAddress("CaloClusters.hits_begin", CaloClusters_hits_begin, &b_CaloClusters_hits_begin);
   fChain->SetBranchAddress("CaloClusters.hits_end", CaloClusters_hits_end, &b_CaloClusters_hits_end);
   fChain->SetBranchAddress("CaloClusters#0", &CaloClusters0_, &b_CaloClusters0_);
   fChain->SetBranchAddress("CaloClusters#0.index", &CaloClusters0_index, &b_CaloClusters0_index);
   fChain->SetBranchAddress("CaloClusters#0.collectionID", &CaloClusters0_collectionID, &b_CaloClusters0_collectionID);
   fChain->SetBranchAddress("ECalBarrelCells", &ECalBarrelCells_, &b_ECalBarrelCells_);
   fChain->SetBranchAddress("ECalBarrelCells.core.cellId", ECalBarrelCells_core_cellId, &b_ECalBarrelCells_core_cellId);
   fChain->SetBranchAddress("ECalBarrelCells.core.energy", ECalBarrelCells_core_energy, &b_ECalBarrelCells_core_energy);
   fChain->SetBranchAddress("ECalBarrelCells.core.time", ECalBarrelCells_core_time, &b_ECalBarrelCells_core_time);
   fChain->SetBranchAddress("ECalBarrelCells.core.bits", ECalBarrelCells_core_bits, &b_ECalBarrelCells_core_bits);
   fChain->SetBranchAddress("GenParticles", &GenParticles_, &b_GenParticles_);
   fChain->SetBranchAddress("GenParticles.core.pdgId", GenParticles_core_pdgId, &b_GenParticles_core_pdgId);
   fChain->SetBranchAddress("GenParticles.core.charge", GenParticles_core_charge, &b_GenParticles_core_charge);
   fChain->SetBranchAddress("GenParticles.core.status", GenParticles_core_status, &b_GenParticles_core_status);
   fChain->SetBranchAddress("GenParticles.core.vertex.x", GenParticles_core_vertex_x, &b_GenParticles_core_vertex_x);
   fChain->SetBranchAddress("GenParticles.core.vertex.y", GenParticles_core_vertex_y, &b_GenParticles_core_vertex_y);
   fChain->SetBranchAddress("GenParticles.core.vertex.z", GenParticles_core_vertex_z, &b_GenParticles_core_vertex_z);
   fChain->SetBranchAddress("GenParticles.core.p4.mass", GenParticles_core_p4_mass, &b_GenParticles_core_p4_mass);
   fChain->SetBranchAddress("GenParticles.core.p4.px", GenParticles_core_p4_px, &b_GenParticles_core_p4_px);
   fChain->SetBranchAddress("GenParticles.core.p4.py", GenParticles_core_p4_py, &b_GenParticles_core_p4_py);
   fChain->SetBranchAddress("GenParticles.core.p4.pz", GenParticles_core_p4_pz, &b_GenParticles_core_p4_pz);
   fChain->SetBranchAddress("GenParticles.core.bits", GenParticles_core_bits, &b_GenParticles_core_bits);
   fChain->SetBranchAddress("GenParticles#0", &GenParticles0_, &b_GenParticles0_);
   fChain->SetBranchAddress("GenParticles#0.index", GenParticles0_index, &b_GenParticles0_index);
   fChain->SetBranchAddress("GenParticles#0.collectionID", GenParticles0_collectionID, &b_GenParticles0_collectionID);
   fChain->SetBranchAddress("GenParticles#1", &GenParticles1_, &b_GenParticles1_);
   fChain->SetBranchAddress("GenParticles#1.index", GenParticles1_index, &b_GenParticles1_index);
   fChain->SetBranchAddress("GenParticles#1.collectionID", GenParticles1_collectionID, &b_GenParticles1_collectionID);
   Notify();
}

Bool_t FourteenLayers40GeV::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FourteenLayers40GeV::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FourteenLayers40GeV::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FourteenLayers40GeV_cxx
