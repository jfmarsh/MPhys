#define SixteenLayersAnalysis_cxx
#include "SixteenLayersAnalysis.h"
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

void SixteenLayersAnalysis::Loop()
{   
   // CALO ANALYSIS 
   if (fChain == 0) return;

   // returns fentries 
   Long64_t nentries = fChain->GetEntriesFast();

   // write histograms to root file 
   TFile * output_file = new TFile("CaloHistograms.root", "recreate");
   TTree * root_tree = new TTree("tree", "root-tuple description");
   root_tree->SetDirectory(output_file);

   // HISTOGRAMS           
   // INVARAINT MASS - e+e- 
   TH1F *InvariantMass = new TH1F("invariantMass","Reconstructed Invariant Mass with 2 clusters (16 Layers); Mass [GeV]; Events",150,0,115);
   TH1F *massHisto = new TH1F("massHisto","CaloCluster Energies Sum (16 Layers); Mass [GeV]; Events",150,0,115);        
   gStyle->SetOptStat("ne");                                             
   // PSEUDORAPIDITY    
   TH1F *PseudoWith = new TH1F("PseudoWith","e+e- Pseudorapidity (with clusters); Pseudorapidity; Events", 100,-5,5);           
   TH1F *PseudoWithout = new TH1F("PseudoWithout","e+e- Pseudorapidity (without clusters); Pseudorapidity; Events", 100,-5,5);   
   TH1F *PseudoClusters = new TH1F("PseudoWith","Pseudorapidity of the clusters; Pseudorapidity; Events", 100,-1,1); 
   // cos theta
   TH1F *costheta = new TH1F("cos_theta","Cos theta 12 (16 Layers); Cos theta; ",100,-1,1);  
   // 3 cluster Invariant Mass
   TH1F *InvariantMass3 = new TH1F("invariantMass3","Reconstructed Invariant Mass [3rd cluster > 10 GeV] (16 Layers); Mass [GeV]; Events",150,0,115);
   
   // determining #particles and #events with/without clusters 
   Int_t EntriesWithout = 0;
   Int_t EntriesWith = 0;
   Int_t Events = 0;

   ///////////// LOOP OVER ALL EVENTS //////////////////
   double nEvents = nentries - 2;
   for (Long64_t jentry=0; jentry<nEvents; jentry++) {

      // print event number, #GenParticles and #CaloClusters 
      GetEntry(jentry);
      std::cout << "Event: " << jentry << "\n";
      std::cout << "#GenParticles: " << GenParticles_ << "\n";
      std::cout << "#CaloClusters: " << CaloClusters_ << "\n";

      // positions of top three most energetic clusters 
      double r1 = sqrt(pow(CaloClusters_core_position_x[0],2) + pow(CaloClusters_core_position_y[0],2) + pow(CaloClusters_core_position_z[0],2));
      double r2 = sqrt(pow(CaloClusters_core_position_x[1],2) + pow(CaloClusters_core_position_y[1],2) + pow(CaloClusters_core_position_z[1],2));
      double r3 = sqrt(pow(CaloClusters_core_position_x[2],2) + pow(CaloClusters_core_position_y[2],2) + pow(CaloClusters_core_position_z[2],2));

      // cos theta 12 
      double numerator12 = CaloClusters_core_position_x[0]*CaloClusters_core_position_x[1] + CaloClusters_core_position_y[0]*CaloClusters_core_position_y[1] + CaloClusters_core_position_z[0]*CaloClusters_core_position_z[1];
      double cosTheta12 = numerator12/(r1*r2);
      costheta->Fill(cosTheta12);

      // cos theta 13
      double numerator13 = CaloClusters_core_position_x[0]*CaloClusters_core_position_x[2] + CaloClusters_core_position_y[0]*CaloClusters_core_position_y[2] + CaloClusters_core_position_z[0]*CaloClusters_core_position_z[2];
      double cosTheta13 = numerator13/(r1*r3);

      // cos theta 23
      double numerator23 = CaloClusters_core_position_x[1]*CaloClusters_core_position_x[2] + CaloClusters_core_position_y[1]*CaloClusters_core_position_y[2] + CaloClusters_core_position_z[1]*CaloClusters_core_position_z[2];
      double cosTheta23 = numerator23/(r2*r3);

      //////////////// 2 CLUSTERS /////////////////////////
      // invariant mass 
      double invariantMass = sqrt(2 * CaloClusters_core_energy[0] * CaloClusters_core_energy[1] * (1 - cosTheta12));
      std::cout << "Invariant Mass: " << invariantMass << "\n";
      if(CaloClusters_ > 1){
         InvariantMass->Fill(invariantMass);
      }
      // sum of cluster energies
      double TotalClusterEnergy = CaloClusters_core_energy[0] + CaloClusters_core_energy[1];
      std::cout << "Sum of Energies: " << TotalClusterEnergy << "\n";
      if(CaloClusters_ > 1){
         massHisto->Fill(TotalClusterEnergy);
      }

      //////////////// 3 CLUSTERS /////////////////////////
      double E1E2 = 2 * CaloClusters_core_energy[0] * CaloClusters_core_energy[1] * (1 - cosTheta12);
      double E1E3 = 2 * CaloClusters_core_energy[0] * CaloClusters_core_energy[2] * (1 - cosTheta13);
      double E2E3 = 2 * CaloClusters_core_energy[1] * CaloClusters_core_energy[2] * (1 - cosTheta23);
      double invariantMass3 = sqrt(E1E2 + E1E3 + E2E3);
      // 2 clusters 
      if(CaloClusters_ == 2){
         InvariantMass3->Fill(invariantMass);
      } 
      // if 3 or more clusters 
      if(CaloClusters_ > 2){
         //InvariantMass3->Fill(invariantMass3);
         // set lower bound on the third cluster 
         if(CaloClusters_core_energy[2] > 10){
            InvariantMass3->Fill(invariantMass3);
         }
         if(CaloClusters_core_energy[2] < 10){
            InvariantMass3->Fill(invariantMass);
         }
      } 

      ////////////////// LOOP OVER CALO-CLUSTERS ///////////////
      for (Int_t j = 0; j < CaloClusters_; j++) 
      {// loop over CaloClusters

         // energy of cluster
         double energy = (CaloClusters_core_energy[j]);
         std::cout << j << ": CaloCluster Energy: " << energy << "\n";
         // pseudorapidity of the clusters 
         double psuedo = asinh(CaloClusters_core_position_z[j]/sqrt(pow(CaloClusters_core_position_x[j],2)+pow(CaloClusters_core_position_y[j],2)));
         PseudoClusters->Fill(psuedo);
      } // Loop over CaloClusters 

      // print particles generated in each event
      std::cout << "Particles Generated: " << "\n";

      if(CaloClusters_==0){
         Events = Events +1;
      }

      ////////////////////// LOOP OVER GENPARTICLES /////////////////////////////
      for (Int_t j = 0; j < GenParticles_; j++) 
      {// Loop over GenParticles 

         // total four momentum 
         double momentum = sqrt(pow(GenParticles_core_p4_px[j],2) + pow(GenParticles_core_p4_py[j],2) +pow(GenParticles_core_p4_pz[j],2));
         // tan theta 
         double tantheta = (sqrt(pow(GenParticles_core_p4_px[j],2)+pow(GenParticles_core_p4_py[j],2)))/GenParticles_core_p4_pz[j];
         // define theta for negative tan theta - will redefine pseudorapidity 
         double theta1 = atan(-1.*tantheta);
         // define theta for positive tan theta 
         double theta2 = atan(tantheta);
         // pT
         double TransverseMom = sqrt(pow(GenParticles_core_p4_py[j],2)+pow(GenParticles_core_p4_px[j],2));

         //////// NO CLUSTERS ////////
         if(CaloClusters_ == 0)
         {// loop over no clusters

            EntriesWithout = EntriesWithout + 1;
            // particles generated and their momentum
            std::cout << "PDGID: " << GenParticles_core_pdgId[j] << ", Energy: " << momentum << "\n"; 

            if(GenParticles_core_pdgId[j]==11 || GenParticles_core_pdgId[j]==-11)
            {// loop over electrons and positrons 

               // pseudorapidity 
               if(tantheta<0){
                  PseudoWithout->Fill(log(tan(theta1/2)));
                  }
               else{
                  PseudoWithout->Fill(-1*log(tan(theta2/2)));
               }

            } // Loop over electrons and positrons
         } // Loop over no clusters

         ////////// WITH CLUSTERS ///////////
         else{
            EntriesWith = EntriesWith + 1;
            // particles generated and their momentum 
            std::cout << "PDGID: " << GenParticles_core_pdgId[j] << ", Energy: " << momentum << "\n"; 
            
            if(GenParticles_core_pdgId[j]==11 || GenParticles_core_pdgId[j]==-11)
            {// Loop over electrons and positrons 

               // pseudorapidity 
               if(tantheta<0){
                  PseudoWith->Fill(log(tan(theta1/2)));
               }
               else{
                  PseudoWith->Fill(-1*log(tan(theta2/2)));
               }
            } // Loop over electrons and positrons 
         } // Loop over clusters    
      }// Loop over GenParticles 
      std::cout << "\n" << endl;
   } // Loop over all events
    
   // #events and #entries with/without clusters 
   std::cout << "#Events without Clusters: " << Events << "\n";
   std::cout << "#Entries without Clusters: " << EntriesWithout << "\n";
   std::cout << "#Entries with Clusters: " << EntriesWith << "\n\n";

   //////////// PLOTTING HISTORAMS ////////////////

   // INVARIANT MASS
   TCanvas * c1 = new TCanvas(); InvariantMass->Draw();
   gStyle->SetStatX(0.5);
   // Fit peak 
   TF1 *fit2 = new TF1("fit2","gaus",85,100);
   gStyle->SetOptFit(1111);
   InvariantMass->Fit(fit2,"LRI");
   Double_t mean2 = fit2->GetParameter(1);
   Double_t mean_error2 = fit2->GetParError(1);
   Double_t sigma2 = fit2->GetParameter(2);
   Double_t sigma_error2 = fit2->GetParError(2);
   // print fit results 
   std::cout << "mean: " << mean2 << " +/- " << mean_error2 << "\n"; 
   std::cout << "sigma: " << sigma2 << " +/- " << sigma_error2 << "\n"; 
   // write to text file 
   ofstream NewerFile("InvariantMassFit_16Layers.txt");
   NewerFile << mean2 << " " << mean_error2 << " " << sigma2 << " " << sigma_error2 << "\n\n";
   NewerFile.close();
   InvariantMass->SetDirectory(0);
   c1->SaveAs("InvariantMassFit_16Layers.pdf");

   // INDIVIDUAL CLUSTER SUM OF ENERGIES 
   TCanvas * c2 = new TCanvas(); massHisto->Draw();
   // FITTING GAUSSIAN PEAK 
   TF1 *fit = new TF1("fit","gaus",85,100);
   massHisto->Fit(fit,"LRI");
   Double_t mean = fit->GetParameter(1);
   Double_t mean_error = fit->GetParError(1);
   Double_t sigma = fit->GetParameter(2);
   Double_t sigma_error = fit->GetParError(2);
   // write to text file 
   ofstream MyFile("ClusterEnergySum_16Layers.txt");
   // write to file in format: mean, mean_error, sigma, sigma_error
   MyFile << mean << " " << mean_error << " " << sigma << " " << sigma_error << "\n\n";
   MyFile.close();
   massHisto->SetDirectory(0);
   c2->SaveAs("ClusterEnergySum_16Layers.pdf");

   // PSEUDORAPIDITY
   TCanvas * c3 = new TCanvas(); PseudoWith->Draw();
   c3->SaveAs("ReconAnalysis_PseudoWithClusters.pdf");
   TCanvas * c4 = new TCanvas(); PseudoWithout->Draw();
   c4->SaveAs("ReconAnalysis_PseudoWithoutClusters.pdf");
   PseudoWith->SetDirectory(0); PseudoWithout->SetDirectory(0);
   TCanvas * c5 = new TCanvas(); PseudoClusters->Draw();
   c5->SaveAs("ReconAnalysis_PseudoClusters.pdf");
   PseudoClusters->SetDirectory(0); 

   // CLUSTER ENERGIES 
   TCanvas * c6 = new TCanvas(); costheta->Draw();
   c6->SaveAs("cosTheta_16layers.pdf");
   costheta->SetDirectory(0);

   // INVARIANT MASS - 3 clusters 
   TCanvas * c7 = new TCanvas(); InvariantMass3->Draw();
   // Fit peak 
   TF1 *fit3 = new TF1("fit3","gaus",85,100);
   InvariantMass3->Fit(fit3,"LRI");
   Double_t mean3 = fit3->GetParameter(1);
   Double_t mean_error3 = fit3->GetParError(1);
   Double_t sigma3 = fit3->GetParameter(2);
   Double_t sigma_error3 = fit3->GetParError(2);
   // write to text file 
   ofstream NewestFile("InvariantMassFit_16Layers_3ClustersBound.txt");
   NewestFile << mean3 << " " << mean_error3 << " " << sigma3 << " " << sigma_error3 << "\n\n";
   NewestFile.close();
   c7->SaveAs("InvariantMassFit_16Layers_3ClustersBound.pdf");
   InvariantMass3->SetDirectory(0);

   // write to output root file 
   output_file->cd();
   massHisto->Write();
   PseudoWith->Write(); PseudoWithout->Write();
   output_file->Close();

}
