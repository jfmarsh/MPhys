#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector> 
using namespace std;
// root stuff
#include <TROOT.h>
#include <TChain.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TFile.h>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TH2.h>

int noiseValueHistograms(void)
{
    // create output file to hold histograms 
    TFile * output_file = new TFile("noise25layers.root", "recreate");
    TTree * root_tree = new TTree("tree", "root-tuple description");
    root_tree->SetDirectory(output_file);

    // set ymin and ymax for histogram range
    double ymax = 0.045;
    double ymin = 0.0;
  
    // create histogram for each layer 
    TH1 * noiseHist_layer1 = new TH1F("h_elecNoise_fcc_1", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer1->SetMaximum(ymax); noiseHist_layer1->SetMinimum(ymin);
    TH1 * noiseHist_layer2 = new TH1F("h_elecNoise_fcc_2", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer2->SetMaximum(ymax); noiseHist_layer2->SetMinimum(ymin);
    TH1 * noiseHist_layer3 = new TH1F("h_elecNoise_fcc_3", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer3->SetMaximum(ymax); noiseHist_layer3->SetMinimum(ymin);
    TH1 * noiseHist_layer4 = new TH1F("h_elecNoise_fcc_4", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer4->SetMaximum(ymax); noiseHist_layer4->SetMinimum(ymin);
    TH1 * noiseHist_layer5 = new TH1F("h_elecNoise_fcc_5", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer5->SetMaximum(ymax); noiseHist_layer5->SetMinimum(ymin);
    TH1 * noiseHist_layer6 = new TH1F("h_elecNoise_fcc_6", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer6->SetMaximum(ymax); noiseHist_layer6->SetMinimum(ymin);
    TH1 * noiseHist_layer7 = new TH1F("h_elecNoise_fcc_7", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer7->SetMaximum(ymax); noiseHist_layer7->SetMinimum(ymin);
    TH1 * noiseHist_layer8 = new TH1F("h_elecNoise_fcc_8", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer8->SetMaximum(ymax); noiseHist_layer8->SetMinimum(ymin);
    TH1 * noiseHist_layer9 = new TH1F("h_elecNoise_fcc_9", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer9->SetMaximum(ymax); noiseHist_layer9->SetMinimum(ymin);
    TH1 * noiseHist_layer10 = new TH1F("h_elecNoise_fcc_10", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer10->SetMaximum(ymax); noiseHist_layer10->SetMinimum(ymin);
    TH1 * noiseHist_layer11 = new TH1F("h_elecNoise_fcc_11", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer11->SetMaximum(ymax); noiseHist_layer11->SetMinimum(ymin);
    TH1 * noiseHist_layer12 = new TH1F("h_elecNoise_fcc_12", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer12->SetMaximum(ymax); noiseHist_layer12->SetMinimum(ymin);
    TH1 * noiseHist_layer13 = new TH1F("h_elecNoise_fcc_13", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer13->SetMaximum(ymax); noiseHist_layer13->SetMinimum(ymin);
    TH1 * noiseHist_layer14 = new TH1F("h_elecNoise_fcc_14", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer14->SetMaximum(ymax); noiseHist_layer14->SetMinimum(ymin);
    TH1 * noiseHist_layer15 = new TH1F("h_elecNoise_fcc_15", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer15->SetMaximum(ymax); noiseHist_layer15->SetMinimum(ymin);
    TH1 * noiseHist_layer16 = new TH1F("h_elecNoise_fcc_16", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer16->SetMaximum(ymax); noiseHist_layer16->SetMinimum(ymin);
    TH1 * noiseHist_layer17 = new TH1F("h_elecNoise_fcc_17", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer17->SetMaximum(ymax); noiseHist_layer17->SetMinimum(ymin);
    TH1 * noiseHist_layer18 = new TH1F("h_elecNoise_fcc_18", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer18->SetMaximum(ymax); noiseHist_layer18->SetMinimum(ymin);
    TH1 * noiseHist_layer19 = new TH1F("h_elecNoise_fcc_19", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer19->SetMaximum(ymax); noiseHist_layer19->SetMinimum(ymin);
    TH1 * noiseHist_layer20 = new TH1F("h_elecNoise_fcc_20", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer20->SetMaximum(ymax); noiseHist_layer20->SetMinimum(ymin);
    TH1 * noiseHist_layer21 = new TH1F("h_elecNoise_fcc_21", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer21->SetMaximum(ymax); noiseHist_layer21->SetMinimum(ymin);
    TH1 * noiseHist_layer22 = new TH1F("h_elecNoise_fcc_22", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer22->SetMaximum(ymax); noiseHist_layer22->SetMinimum(ymin);
    TH1 * noiseHist_layer23 = new TH1F("h_elecNoise_fcc_23", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer23->SetMaximum(ymax); noiseHist_layer23->SetMinimum(ymin);
    TH1 * noiseHist_layer24 = new TH1F("h_elecNoise_fcc_24", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer24->SetMaximum(ymax); noiseHist_layer24->SetMinimum(ymin);
    TH1 * noiseHist_layer25 = new TH1F("h_elecNoise_fcc_25", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer25->SetMaximum(ymax); noiseHist_layer25->SetMinimum(ymin);
    TH1 * noiseHist_layer26 = new TH1F("h_elecNoise_fcc_26", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer26->SetMaximum(ymax); noiseHist_layer26->SetMinimum(ymin);
    TH1 * noiseHist_layer27 = new TH1F("h_elecNoise_fcc_27", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer27->SetMaximum(ymax); noiseHist_layer27->SetMinimum(ymin);
    TH1 * noiseHist_layer28 = new TH1F("h_elecNoise_fcc_28", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer28->SetMaximum(ymax); noiseHist_layer28->SetMinimum(ymin);
    TH1 * noiseHist_layer29 = new TH1F("h_elecNoise_fcc_29", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer29->SetMaximum(ymax); noiseHist_layer29->SetMinimum(ymin);
    TH1 * noiseHist_layer30 = new TH1F("h_elecNoise_fcc_30", "Electronic Noise; |eta|; Electronic Noise [GeV]", 83, 0, 0.825);
    noiseHist_layer30->SetMaximum(ymax); noiseHist_layer30->SetMinimum(ymin);

    /////////// READ IN CONSTANT FIT VALUES ///////////
    // create an array to hold the values 
    vector<double> constant;
    // create an input file stream 
    ifstream file("FlatFitResults.txt",ios::in);
    // variable to hold each number as it reads
    double number;
    // read number using the extraction (>>) operator
    while (file >> number){
        // add the number to the end of the array 
        constant.push_back(number);
    }
    // Close the file 
    file.close();
  
    /////////// READ IN LINEAR FIT VALUES ///////////
    // GRADIENT VALUES 
    vector<double> Gradient;
    ifstream file2("LinearFitResults_Gradient.txt",ios::in);
    while (file2 >> number){
        Gradient.push_back(number);
    }
    file2.close();

    // INTERCEPT VALUES 
    vector<double> Intercept;
    ifstream file3("LinearFitResults_Intercept.txt",ios::in);
    while (file3 >> number){
        Intercept.push_back(number);
    }
    file3.close();

    /////////// ETA VALUES  ///////////
    vector<double> etaValues;
    ifstream file4("etaValues.txt",ios::in);
    while (file4 >> number){
        etaValues.push_back(number);
    }
    file4.close();

    /////////// X VALUES - # of layers //////////
    int num_layers = 25;
    // merging of layer defined as: [a]*1 + [b]*num_layers-2 + [c]*1
    int a = 4;
    int b = 3;
    int c = 17;
    // creating vector containing x values 
    std::vector <int> layers (num_layers);
    std::cout << "Number of layers: " << layers.size() << "\n";
    for(int i=0; i<layers.size(); i++)
    {
        if(i==0){
            layers[i] = 1;
        }
        if(i==1){
            layers[i] = 1 + a;
        }
        if(i==(num_layers-1)){
            layers[i] = 90 - c;
        }
    }
    int interval = layers.size() - 1;
    for(int i=2; i<interval; i++){
        layers[i] = layers[i-1] + b;
        
    }
    for(int i=0; i<layers.size(); i++)
    {
       cout << layers[i] << "\n";
    }

    for (int j=0; j<etaValues.size(); j++)
    { // loop over eta 

        double etaVal = etaValues[j]; 
        
        // Fill histograms with constant noise 
        noiseHist_layer1->SetBinContent(noiseHist_layer1->FindBin(etaVal), constant[j]);
        noiseHist_layer2->SetBinContent(noiseHist_layer2->FindBin(etaVal), constant[j]);
        noiseHist_layer3->SetBinContent(noiseHist_layer3->FindBin(etaVal), constant[j]);
        noiseHist_layer4->SetBinContent(noiseHist_layer4->FindBin(etaVal), constant[j]);
        noiseHist_layer5->SetBinContent(noiseHist_layer5->FindBin(etaVal), constant[j]);
        noiseHist_layer6->SetBinContent(noiseHist_layer6->FindBin(etaVal), constant[j]);
        noiseHist_layer7->SetBinContent(noiseHist_layer7->FindBin(etaVal), constant[j]);
        noiseHist_layer8->SetBinContent(noiseHist_layer8->FindBin(etaVal), constant[j]);
        noiseHist_layer9->SetBinContent(noiseHist_layer9->FindBin(etaVal), constant[j]);
        
        // LINEAR function - noise = gradient * #layers + intercept 
        double noise10 = Gradient[j] * layers[9] + Intercept[j]; noiseHist_layer10->SetBinContent(noiseHist_layer10->FindBin(etaVal), noise10);
        double noise11 = Gradient[j] * layers[10] + Intercept[j]; noiseHist_layer11->SetBinContent(noiseHist_layer11->FindBin(etaVal), noise11);
        double noise12 = Gradient[j] * layers[11] + Intercept[j]; noiseHist_layer12->SetBinContent(noiseHist_layer12->FindBin(etaVal), noise12);
        double noise13 = Gradient[j] * layers[12] + Intercept[j]; noiseHist_layer13->SetBinContent(noiseHist_layer13->FindBin(etaVal), noise13);
        double noise14 = Gradient[j] * layers[13] + Intercept[j]; noiseHist_layer14->SetBinContent(noiseHist_layer14->FindBin(etaVal), noise14);
        double noise15 = Gradient[j] * layers[14] + Intercept[j]; noiseHist_layer15->SetBinContent(noiseHist_layer15->FindBin(etaVal), noise15);
        double noise16 = Gradient[j] * layers[15] + Intercept[j]; noiseHist_layer16->SetBinContent(noiseHist_layer16->FindBin(etaVal), noise16);
        double noise17 = Gradient[j] * layers[16] + Intercept[j]; noiseHist_layer17->SetBinContent(noiseHist_layer17->FindBin(etaVal), noise17);
        double noise18 = Gradient[j] * layers[17] + Intercept[j]; noiseHist_layer18->SetBinContent(noiseHist_layer18->FindBin(etaVal), noise18);
        double noise19 = Gradient[j] * layers[18] + Intercept[j]; noiseHist_layer19->SetBinContent(noiseHist_layer19->FindBin(etaVal), noise19);
        double noise20 = Gradient[j] * layers[19] + Intercept[j]; noiseHist_layer20->SetBinContent(noiseHist_layer20->FindBin(etaVal), noise20);
        double noise21 = Gradient[j] * layers[20] + Intercept[j]; noiseHist_layer21->SetBinContent(noiseHist_layer21->FindBin(etaVal), noise21);
        double noise22 = Gradient[j] * layers[21] + Intercept[j]; noiseHist_layer22->SetBinContent(noiseHist_layer22->FindBin(etaVal), noise22);
        double noise23 = Gradient[j] * layers[22] + Intercept[j]; noiseHist_layer23->SetBinContent(noiseHist_layer23->FindBin(etaVal), noise23);
        double noise24 = Gradient[j] * layers[23] + Intercept[j]; noiseHist_layer24->SetBinContent(noiseHist_layer24->FindBin(etaVal), noise24);
        double noise25 = Gradient[j] * layers[24] + Intercept[j]; noiseHist_layer25->SetBinContent(noiseHist_layer25->FindBin(etaVal), noise25);
        double noise26 = Gradient[j] * layers[25] + Intercept[j]; noiseHist_layer26->SetBinContent(noiseHist_layer26->FindBin(etaVal), noise26);
        double noise27 = Gradient[j] * layers[26] + Intercept[j]; noiseHist_layer27->SetBinContent(noiseHist_layer27->FindBin(etaVal), noise27);
        double noise28 = Gradient[j] * layers[27] + Intercept[j]; noiseHist_layer28->SetBinContent(noiseHist_layer28->FindBin(etaVal), noise28);
        double noise29 = Gradient[j] * layers[28] + Intercept[j]; noiseHist_layer29->SetBinContent(noiseHist_layer29->FindBin(etaVal), noise29);
        double noise30 = Gradient[j] * layers[29] + Intercept[j]; noiseHist_layer30->SetBinContent(noiseHist_layer30->FindBin(etaVal), noise30);

    } // loop over eta 
  
    // Draw each histograms on a new canvas 
    TCanvas * c1 = new TCanvas(); noiseHist_layer1->Draw();
    TCanvas * c2 = new TCanvas(); noiseHist_layer2->Draw();
    TCanvas * c3 = new TCanvas(); noiseHist_layer3->Draw();
    TCanvas * c4 = new TCanvas(); noiseHist_layer4->Draw();
    TCanvas * c5 = new TCanvas(); noiseHist_layer5->Draw();
    TCanvas * c6 = new TCanvas(); noiseHist_layer6->Draw();
    TCanvas * c7 = new TCanvas(); noiseHist_layer7->Draw();
    TCanvas * c8 = new TCanvas(); noiseHist_layer8->Draw();
    TCanvas * c9 = new TCanvas(); noiseHist_layer9->Draw();
    TCanvas * c10 = new TCanvas(); noiseHist_layer10->Draw();
    TCanvas * c11 = new TCanvas(); noiseHist_layer11->Draw();
    TCanvas * c12 = new TCanvas(); noiseHist_layer12->Draw();
    TCanvas * c13 = new TCanvas(); noiseHist_layer13->Draw();
    TCanvas * c14 = new TCanvas(); noiseHist_layer14->Draw();
    TCanvas * c15 = new TCanvas(); noiseHist_layer15->Draw();
    TCanvas * c16 = new TCanvas(); noiseHist_layer16->Draw();
    TCanvas * c17 = new TCanvas(); noiseHist_layer17->Draw();
    TCanvas * c18 = new TCanvas(); noiseHist_layer18->Draw();
    TCanvas * c19 = new TCanvas(); noiseHist_layer19->Draw();
    TCanvas * c20 = new TCanvas(); noiseHist_layer20->Draw();
    TCanvas * c21 = new TCanvas(); noiseHist_layer21->Draw();
    TCanvas * c22 = new TCanvas(); noiseHist_layer22->Draw();
    TCanvas * c23 = new TCanvas(); noiseHist_layer23->Draw();
    TCanvas * c24 = new TCanvas(); noiseHist_layer24->Draw();
    TCanvas * c25 = new TCanvas(); noiseHist_layer25->Draw();
    TCanvas * c26 = new TCanvas(); noiseHist_layer26->Draw();
    TCanvas * c27 = new TCanvas(); noiseHist_layer27->Draw();
    TCanvas * c28 = new TCanvas(); noiseHist_layer28->Draw();
    TCanvas * c29 = new TCanvas(); noiseHist_layer29->Draw();
    TCanvas * c30 = new TCanvas(); noiseHist_layer30->Draw();
  
    noiseHist_layer1->SetDirectory(0); noiseHist_layer2->SetDirectory(0);
    noiseHist_layer3->SetDirectory(0); noiseHist_layer4->SetDirectory(0);
    noiseHist_layer5->SetDirectory(0); noiseHist_layer6->SetDirectory(0);
    noiseHist_layer7->SetDirectory(0); noiseHist_layer8->SetDirectory(0);
    noiseHist_layer9->SetDirectory(0); noiseHist_layer10->SetDirectory(0);
    noiseHist_layer11->SetDirectory(0); noiseHist_layer12->SetDirectory(0);
    noiseHist_layer13->SetDirectory(0); noiseHist_layer14->SetDirectory(0);
    noiseHist_layer15->SetDirectory(0); noiseHist_layer16->SetDirectory(0);
    noiseHist_layer17->SetDirectory(0); noiseHist_layer18->SetDirectory(0);
    noiseHist_layer19->SetDirectory(0); noiseHist_layer20->SetDirectory(0);
    noiseHist_layer21->SetDirectory(0); noiseHist_layer22->SetDirectory(0);
    noiseHist_layer23->SetDirectory(0); noiseHist_layer24->SetDirectory(0);
    noiseHist_layer25->SetDirectory(0); noiseHist_layer26->SetDirectory(0);
    noiseHist_layer27->SetDirectory(0); noiseHist_layer28->SetDirectory(0);
    noiseHist_layer29->SetDirectory(0); noiseHist_layer30->SetDirectory(0);

    // write histograms to output file 
    output_file->cd();
    noiseHist_layer1->Write(); noiseHist_layer2->Write();
    noiseHist_layer3->Write(); noiseHist_layer4->Write();
    noiseHist_layer5->Write(); noiseHist_layer6->Write();
    noiseHist_layer7->Write(); noiseHist_layer8->Write();
    noiseHist_layer9->Write(); noiseHist_layer10->Write();
    noiseHist_layer11->Write(); noiseHist_layer12->Write();
    noiseHist_layer13->Write(); noiseHist_layer14->Write();
    noiseHist_layer15->Write(); noiseHist_layer16->Write();
    noiseHist_layer17->Write(); noiseHist_layer18->Write();
    noiseHist_layer19->Write(); noiseHist_layer20->Write();
    noiseHist_layer21->Write(); noiseHist_layer22->Write();
    noiseHist_layer23->Write(); noiseHist_layer24->Write();
    noiseHist_layer25->Write(); noiseHist_layer26->Write();
    noiseHist_layer27->Write(); noiseHist_layer28->Write();
    noiseHist_layer29->Write(); noiseHist_layer30->Write();

    output_file->Close();
    
  return 0;

}