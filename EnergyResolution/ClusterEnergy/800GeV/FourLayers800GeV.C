#define FourLayers800GeV_cxx
#include "FourLayers800GeV.h"
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

void FourLayers800GeV::Loop()
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
   TH1F *ClusterEnergy = new TH1F("clusterEnergy"," CaloCluster Energy for 800 GeV Electron Beam (4 Layers); Energy [GeV]; #Events",150,400,850);    
   gStyle->SetOptStat("ne");   
   TH1F *EnergyPeak = new TH1F("EnergyPeak"," CaloCluster Energy for 800 GeV Electron Beam (4 Layers); Energy [GeV]; #Events",150,400,850);                                 

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
            if(CaloClusters_core_energy[0] < 798){
               totalEnergy = CaloClusters_core_energy[0] + CaloClusters_core_energy[1];
            }
            if(CaloClusters_core_energy[0] > 798){
               totalEnergy = CaloClusters_core_energy[0];
            }
         }
      } // loop over caloclutsers 
      ClusterEnergy->Fill(totalEnergy);
      //ClusterEnergy->Fill(CaloClusters_core_energy[0]);
      EnergyPeak->Fill(CaloClusters_core_energy[0]);
   } // Loop over all events
    

   //////////// PLOTTING HISTORAMS ////////////////

   double maxVal = ClusterEnergy->GetXaxis()->GetBinCenter(ClusterEnergy->GetMaximumBin());
   std::cout << "\nMax Value: " << maxVal << "\n\n" << endl;

   //  Cluster energy 
   TCanvas * c1 = new TCanvas(); ClusterEnergy->Draw();
   gStyle->SetStatX(0.5);
   // Fit peak 
   TF1 *fit = new TF1("fit","gaus",740,840);
   gStyle->SetOptFit(1011);
   ClusterEnergy->Fit(fit,"LRI");
   Double_t mean = fit->GetParameter(1);
   Double_t mean_error = fit->GetParError(1);
   Double_t sigma = fit->GetParameter(2);
   Double_t sigma_error = fit->GetParError(2);
   // write to text file 
   ofstream NewerFile("ClusterEnergyFit_800GeV_4Layers_BOUND.txt");
   NewerFile << mean << " " << mean_error << " " << sigma << " " << sigma_error << "\n\n";
   NewerFile.close();
   c1->SaveAs("ClusterEnergyFit_800GeV_4layers_BOUND.pdf");
   ClusterEnergy->SetDirectory(0);

   //  TWO GAUSSIANS
   TCanvas * c2 = new TCanvas(); EnergyPeak->Draw();
   gStyle->SetStatX(0.5);
   // create fits 
   TF1 *fit1 = new TF1("fit1","gaus",730, maxVal);
   TF1 *fit2 = new TF1("fit2","gaus",maxVal,840);
   TF1 *total = new TF1("total", "gaus(0)+gaus(3)",730, 825);
   gStyle->SetOptFit(1011);
   // parameter array for the total function 
   double par[6];
   fit1->SetLineColor(4); 
   EnergyPeak->Fit(fit1,"LIR+");
   EnergyPeak->Fit(fit2,"LIR+");
   // get the parameters from the fit 
   fit1->GetParameters(&par[0]);
   fit2->GetParameters(&par[3]);
   total->SetParameters(par);
   //fit1->SetLineColor(4); 
   //EnergyPeak->Fit("total","LIR+");

   // extract fit parameters 
   Double_t mean1 = fit1->GetParameter(1);
   Double_t mean_error1 = fit1->GetParError(1);
   Double_t sigma1 = fit1->GetParameter(2);
   Double_t sigma_error1 = fit1->GetParError(2);
   // write to text file 
   ofstream File1("ClusterEnergyFit_800GeV_4Layers_gauss1.txt");
   File1 << mean1 << " " << mean_error1 << " " << sigma1 << " " << sigma_error1 << "\n\n";
   File1.close();

   Double_t mean2 = fit2->GetParameter(1);
   Double_t mean_error2 = fit2->GetParError(1);
   Double_t sigma2 = fit2->GetParameter(2);
   Double_t sigma_error2 = fit2->GetParError(2);
   // write to text file 
   ofstream File2("ClusterEnergyFit_800GeV_4Layers_gauss2.txt");
   File2 << mean2 << " " << mean_error2 << " " << sigma2 << " " << sigma_error2 << "\n\n";
   File2.close();

   ofstream File3("ClusterEnergyFit_800GeV_4Layers_totalGauss.txt");
   File3 << total->GetParameter(0) << " " << total->GetParameter(1) << " " << total->GetParameter(2) << " " << total->GetParameter(3) << " " << \
   total->GetParameter(4) << " " << total->GetParameter(5) <<  endl;
   //File3 << mean3 << " " << mean_error3 << " " << sigma3 << " " << sigma_error3 << "\n\n";
   File3.close();
   
   c2->SaveAs("ClusterEnergyFit_800GeV_4layers_2Gauss.pdf");
   EnergyPeak->SetDirectory(0);

   // write to output root file 
   output_file->cd();
   ClusterEnergy->Write();
   EnergyPeak->Write();
   output_file->Close();

}
