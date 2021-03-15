#define FourteenLayers400GeV_cxx
#include "FourteenLayers400GeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>

void FourteenLayers400GeV::Loop()
{  // CALO ANALYSIS 
   if (fChain == 0) return;

   // returns fentries 
   Long64_t nentries = fChain->GetEntriesFast();

   // write histograms to root file 
   TFile * output_file = new TFile("CaloHistograms.root", "recreate");
   TTree * root_tree = new TTree("tree", "root-tuple description");
   root_tree->SetDirectory(output_file);

   // HISTOGRAMS           
   // Cluster Energy
   TH1F *ClusterEnergy = new TH1F("clusterEnergy"," CaloCluster Energy for 400 GeV Electron Beam (14 Layers); Energy [GeV]; #Events",150,0,450);    
   gStyle->SetOptStat("ne");                                

   ///////////// LOOP OVER ALL EVENTS //////////////////
   int nEvents = nentries - 2;
   for (Long64_t jentry=0; jentry<nEvents; jentry++) {

      // print event number, #GenParticles and #CaloClusters 
      GetEntry(jentry);
      std::cout << "Event: " << jentry << "\n";
      std::cout << "#GenParticles: " << GenParticles_ << "\n";
      std::cout << "#CaloClusters: " << CaloClusters_ << "\n";

      double totalEnergy = 0.0; 
      for (Int_t j = 0; j < CaloClusters_; j++) 
      {// loop over CaloClusters
         double energy = (CaloClusters_core_energy[j]);
         std::cout << j << ": CaloCluster Energy: " << energy << "\n";
      } // loop over caloclutsers 
      ClusterEnergy->Fill(CaloClusters_core_energy[0]);
   } // Loop over all events

   //////////// PLOTTING HISTORAMS ////////////////

   //  Cluster energy 
   TCanvas * c1 = new TCanvas(); ClusterEnergy->Draw();
   gStyle->SetStatX(0.5);
   // Fit peak 
   TF1 *fit = new TF1("fit","gaus",385,410);
   gStyle->SetOptFit(1011);
   ClusterEnergy->Fit(fit,"LRI");
   Double_t mean = fit->GetParameter(1);
   Double_t mean_error = fit->GetParError(1);
   Double_t sigma = fit->GetParameter(2);
   Double_t sigma_error = fit->GetParError(2);
   // write to text file 
   ofstream NewerFile("ClusterEnergyFit_400GeV_14Layers.txt");
   NewerFile << mean << " " << mean_error << " " << sigma << " " << sigma_error << "\n\n";
   NewerFile.close();
   c1->SaveAs("ClusterEnergyFit_400GeV_14layers.pdf");
   ClusterEnergy->SetDirectory(0);

   // write to output root file 
   output_file->cd();
   ClusterEnergy->Write();
   output_file->Close();

}

