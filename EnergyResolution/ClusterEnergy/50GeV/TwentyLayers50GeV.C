#define TwentyLayers50GeV_cxx
#include "TwentyLayers50GeV.h"
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

void TwentyLayers50GeV::Loop()
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
   TH1F *ClusterEnergy = new TH1F("clusterEnergy"," CaloCluster Energy for 50 GeV Electron Beam (20 Layers); Energy [GeV]; #Clusters",150,0,60);    
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

         if(CaloClusters_ == 1){
            totalEnergy = CaloClusters_core_energy[j];
         }
         if(CaloClusters_ > 1){
            if(CaloClusters_core_energy[0] < 49){
               totalEnergy = CaloClusters_core_energy[0] + CaloClusters_core_energy[1];
            }
            if(CaloClusters_core_energy[0] > 49){
               totalEnergy = CaloClusters_core_energy[0];
            }
         }
      } // loop over caloclutsers 
      ClusterEnergy->Fill(totalEnergy);
      //ClusterEnergy->Fill(CaloClusters_core_energy[0]);
   } // Loop over all events
    

   //////////// PLOTTING HISTORAMS ////////////////

   //  Cluster energy 
   TCanvas * c1 = new TCanvas(); ClusterEnergy->Draw();
   gStyle->SetStatX(0.5);
   // Fit peak 
   TF1 *fit = new TF1("fit","gaus",46,54);
   gStyle->SetOptFit(1011);
   ClusterEnergy->Fit(fit,"LRI");
   Double_t mean = fit->GetParameter(1);
   Double_t mean_error = fit->GetParError(1);
   Double_t sigma = fit->GetParameter(2);
   Double_t sigma_error = fit->GetParError(2);
   // write to text file 
   ofstream NewerFile("ClusterEnergyFit_50GeV_20Layers_BOUND.txt");
   NewerFile << mean << " " << mean_error << " " << sigma << " " << sigma_error << "\n\n";
   NewerFile.close();
   c1->SaveAs("ClusterEnergyFit_50GeV_20layers_BOUND.pdf");
   ClusterEnergy->SetDirectory(0);

   // write to output root file 
   output_file->cd();
   ClusterEnergy->Write();
   output_file->Close();

}

