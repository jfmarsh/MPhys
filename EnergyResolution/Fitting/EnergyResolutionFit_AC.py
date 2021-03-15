from __future__ import print_function
import ROOT
import numpy as np
from ROOT import TCanvas, TGraph
from ROOT import TF1
from ROOT import gROOT
from array import array
from ROOT import TLegend
from ROOT import TGraphErrors

# Function to calculate the energy resolution
def energyRes(sigma, mean):
    y = sigma/mean
    return y 

# Function to calculate the error on the energy resolution
def y_error(y, par):
    # par is an array of the form [mean, mean_error, sigma, sigma_error]
    delta_y = y * np.sqrt((par[3]/par[2])**2 + (par[1]/par[0])**2)
    return delta_y

layers = [4, 8, 12, 14, 15, 20]
energy = [25, 40, 50, 100, 300, 350, 400]
# open dictionaries to contain fit parameters and errors
energyRes_layer = {}
errors_layer = {}

# loop over layers
for j in layers:
    energyRes_layer[j] = []
    errors_layer[j] = []
    # loop over energy values 
    for E in energy:
        layer = j
        # open text file 
        results = np.loadtxt(f'ClusterEnergy/ClusterEnergyFit_{E}GeV_{layer}Layers.txt')
        # resolution
        resolution = energyRes(results[2], results[0])
        error = y_error(resolution, results)
        # append values to dictionaries 
        energyRes_layer[j].append(resolution)
        errors_layer[j].append(error)

# loop for 200 GeV - layer 15 did not work 
layers_new = [4, 8, 12, 14, 20]
for i in layers_new:
    layer = i
    E = 200
    # open text file 
    results = np.loadtxt(f'ClusterEnergy/ClusterEnergyFit_{E}GeV_{layer}Layers.txt')
    # resolution
    resolution = energyRes(results[2], results[0])
    error = y_error(resolution, results)
    # append values to dictionaries 
    energyRes_layer[i].append(resolution)
    errors_layer[i].append(error)

# define energy values for the plots 
energies = [25, 40, 50, 100, 300, 350, 400, 200]
energy15 = [25, 40, 50, 100, 300, 350, 400]

# create canvas 
c1 = TCanvas( 'c1', 'Energy Resolution', 500, 10, 800, 500 )

# loop over 5 energy values 
n = 8
x4, y4, e4 = array( 'd' ), array( 'd' ), array('d')
x8, y8, e8 = array( 'd' ), array( 'd' ), array('d')
x12, y12, e12 = array( 'd' ), array( 'd' ), array( 'd' )
x14, y14, e14 = array( 'd' ), array( 'd' ), array( 'd' )
x20, y20, e20  = array( 'd' ), array( 'd' ), array('d')
for i in range(len(energies)):
    x4.append( energies[i])
    y4.append( energyRes_layer[4][i])
    e4.append( errors_layer[4][i])
    x8.append( energies[i])
    y8.append( energyRes_layer[8][i])
    e8.append( errors_layer[8][i])
    x12.append( energies[i])
    y12.append( energyRes_layer[12][i])
    e12.append( errors_layer[12][i])
    x14.append( energies[i])
    y14.append( energyRes_layer[14][i])
    e14.append( errors_layer[14][i])
    x20.append( energies[i])
    y20.append( energyRes_layer[20][i])
    e20.append( errors_layer[20][i])
   
# loop over 4 energy values for 15 layers 
n2 = 7
x15, y15, e15 = array( 'd' ), array( 'd' ), array( 'd' )
for i in range(len(energy15) ):
    x15.append( energy15[i])
    y15.append( energyRes_layer[15][i])
    e15.append( errors_layer[15][i])

# LAYER 4
gr = TGraphErrors( n, x4, y4, 0, e4)
gr.SetMarkerColor( 4 )
gr.SetMarkerStyle( 21 )
gr.SetLineColor(4)
# fit 
f1 = TF1("f1","sqrt([0]**2/x + [1]**2)",0,400)
f1.SetParameter(0,1)
f1.SetParameter(1,0.5)
f1.SetLineColor( 4 )
result = gr.Fit(f1, "SQR")
gr.SetTitle( 'Energy Resolution (a and c only)' )
gr.GetXaxis().SetTitle( 'E [GeV]' )
gr.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr.GetYaxis().SetLimits(0,0.2)
gr.Draw( 'AP' )

# LAYER 8
gr2 = TGraphErrors( n, x8, y8, 0, e8)
gr2.SetMarkerColor( 6 )
gr2.SetMarkerStyle( 31 )
gr2.SetLineColor(6)
# fit
f2 = TF1("f2","sqrt([0]**2/x + [1]**2)",0,400)
f2.SetParameter(0,1)
f2.SetParameter(1,0.5)
f2.SetLineColor( 6 )
result2 = gr2.Fit(f2, "SQR")
gr2.Draw('P')

# LAYER 12
gr3 = TGraphErrors( n, x12, y12, 0, e12)
gr3.SetMarkerColor( 80 )
gr3.SetMarkerStyle( 45 )
gr3.SetLineColor(80)
# fit
#f3 = TF1("f3","[0]/sqrt(x) + [1]/x + [2]",0,400)
f3 = TF1("f3","sqrt([0]**2/x + [1]**2)",0,400)
f3.SetParameter(0,1)
f3.SetParameter(1,0.5)
f3.SetLineColor( 80 )
result3 = gr3.Fit(f3, "SQR")
gr3.Draw('P')

# LAYER 15
gr4 = TGraphErrors( n2, x15, y15, 0, e15)
gr4.SetMarkerColor( 2 )
gr4.SetMarkerStyle( 34 )
gr4.SetLineColor(2)
# fit
f4 = TF1("f4","sqrt([0]**2/x + [1]**2)",0,400)
f4.SetParameter(0,1)
f4.SetParameter(1,0.5)
f4.SetLineColor( 2 )
result4 = gr4.Fit(f4, "SQR")
gr4.Draw('P')

# LAYER 20
gr5 = TGraphErrors( n, x20, y20, 0, e20)
gr5.SetMarkerColor( 65 )
gr5.SetMarkerStyle( 47 )
gr5.SetLineColor(65)
# fit
f5 = TF1("f5","sqrt([0]**2/x + [1]**2)",0,400)
f5.SetParameter(0,1)
f5.SetParameter(1,0.5)
f5.SetLineColor( 65 )
result5 = gr5.Fit(f5, "SQR")
gr5.Draw('P')

# LAYER 14
gr6 = TGraphErrors( n, x14, y14, 0, e14)
gr6.SetMarkerColor( 90 )
gr6.SetMarkerStyle( 8 )
gr6.SetLineColor(90)
# fit
f6 = TF1("f6","sqrt([0]**2/x + [1]**2)",0,400)
f6.SetParameter(0,1)
f6.SetParameter(1,0.5)
f6.SetLineColor( 90 )
result6 = gr6.Fit(f6, "SQR")
gr6.Draw('P')

legend = TLegend(0.9,0.7,0.75,0.9)
legend.SetHeader("Legend","C")
legend.AddEntry(gr,"4 layers","p")
legend.AddEntry(gr2,"8 layers","p")
legend.AddEntry(gr3,"12 layers","p")
legend.AddEntry(gr6,"14 layers","p")
legend.AddEntry(gr4,"15 layers","p")
legend.AddEntry(gr5,"20 layers","p")
legend.Draw()

c1.SaveAs("FitResults/EnergyResolutionFit_AC_25to400GeV.pdf")

# extract fit results into text files 
fitResults_layer = {}
fitErrors_layer = {}
fitResults_layer[4] = []
fitErrors_layer[4] = []
fitResults_layer[8] = []
fitErrors_layer[8] = []
fitResults_layer[12] = []
fitErrors_layer[12] = []
fitResults_layer[14] = []
fitErrors_layer[14] = []
fitResults_layer[15] = []
fitErrors_layer[15] = []
fitResults_layer[20] = []
fitErrors_layer[20] = []

for i in range(2):

    fitResults_layer[4].append(result.Get().Parameter(i) * 100)
    fitErrors_layer[4].append(result.Get().Error(i) * 100)
    
    fitResults_layer[8].append(result2.Get().Parameter(i) * 100)
    fitErrors_layer[8].append(result2.Get().Error(i) * 100)
    
    fitResults_layer[12].append(result3.Get().Parameter(i) * 100)
    fitErrors_layer[12].append(result3.Get().Error(i) * 100)

    fitResults_layer[14].append(result6.Get().Parameter(i) * 100)
    fitErrors_layer[14].append(result6.Get().Error(i) * 100)
    
    fitResults_layer[15].append(result4.Get().Parameter(i) * 100)
    fitErrors_layer[15].append(result4.Get().Error(i) * 100)
    
    fitResults_layer[20].append(result5.Get().Parameter(i) * 100)
    fitErrors_layer[20].append(result5.Get().Error(i) * 100)
    
with open('FitResults/EnergyResolutionFitResults_AC_25to400GeV.txt', 'w') as f:
    print(fitResults_layer, file=f) 
f.close

with open('FitResults/EnergyResolutionFitErrors_AC_25to400GeV.txt', 'w') as f:
    print(fitErrors_layer, file=f) 
f.close

# individual graphs 

# LAYER 4
ROOT.gStyle.SetOptFit(1011)
c2 = TCanvas( 'c2', 'Energy Resolution', 500, 10, 800, 500 )
gr = TGraphErrors( n, x4, y4, 0, e4)
gr.SetMarkerColor( 4 )
gr.SetMarkerStyle( 21 )
gr.SetLineColor(4)
f1 = TF1("f1","sqrt([0]**2/x + [1]**2)",0,400)
f1.SetParameter(0,1)
f1.SetParameter(1,0.5)
f1.SetLineColor( 4 )
result = gr.Fit(f1, "SQR")
gr.SetTitle( 'Energy Resolution (a and c only)' )
gr.GetXaxis().SetTitle( 'E [GeV]' )
gr.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr.Draw( 'AP' )
legend = TLegend(0.3,0.8,0.5,0.9)
legend.SetHeader("Legend","C")
legend.AddEntry(gr,"4 layers","p")
legend.Draw()
c2.Draw()
c2.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_4layers.pdf")

# LAYER 8
c3 = TCanvas( 'c3', 'Energy Resolution', 500, 10, 800, 500 )
gr2 = TGraphErrors( n, x8, y8, 0, e8)
gr2.SetMarkerColor( 6 )
gr2.SetMarkerStyle( 31 )
gr2.SetLineColor(6)
f2 = TF1("f2","sqrt([0]**2/x + [1]**2)",0,400)
f2.SetParameter(0,1)
f2.SetParameter(1,0.5)
f2.SetLineColor( 6 )
result2 = gr2.Fit(f2, "SQR")
gr2.SetTitle( 'Energy Resolution (a and c only)' )
gr2.GetXaxis().SetTitle( 'E [GeV]' )
gr2.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr2.Draw( 'AP' )
legend2 = TLegend(0.3,0.8,0.5,0.9)
legend2.SetHeader("Legend","C")
legend2.AddEntry(gr2,"8 layers","p")
legend2.Draw()
c3.Draw()
c3.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_8layers.pdf")

# LAYER 12
c4 = TCanvas( 'c4', 'Energy Resolution', 500, 10, 800, 500 )
gr3 = TGraphErrors( n, x12, y12, 0, e12)
gr3.SetMarkerColor( 80 )
gr3.SetMarkerStyle( 45 )
gr3.SetLineColor(80)
f3 = TF1("f3","sqrt([0]**2/x + [1]**2)",0,400)
f3.SetParameter(0,1)
f3.SetParameter(1,0.5)
f3.SetLineColor( 80 )
result3 = gr3.Fit(f3, "SQR")
gr3.SetTitle( 'Energy Resolution (a and c only)' )
gr3.GetXaxis().SetTitle( 'E [GeV]' )
gr3.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr3.Draw( 'AP' )
legend3 = TLegend(0.3,0.8,0.5,0.9)
legend3.SetHeader("Legend","C")
legend3.AddEntry(gr3,"12 layers","p")
legend3.Draw()
c4.Draw()
c4.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_12layers.pdf")

# LAYER 14
c5 = TCanvas( 'c5', 'Energy Resolution', 500, 10, 800, 500 )
gr6 = TGraphErrors( n, x14, y14, 0, e14)
gr6.SetMarkerColor( 90 )
gr6.SetMarkerStyle( 8 )
gr6.SetLineColor(90)
f6 = TF1("f6","sqrt([0]**2/x + [1]**2)",0,400)
f6.SetParameter(0,1)
f6.SetParameter(1,0.5)
f6.SetLineColor( 90 )
result6 = gr6.Fit(f6, "SQR")
gr6.SetTitle( 'Energy Resolution (a and c only)' )
gr6.GetXaxis().SetTitle( 'E [GeV]' )
gr6.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr6.Draw( 'AP' )
legend6 = TLegend(0.3,0.8,0.5,0.9)
legend6.SetHeader("Legend","C")
legend6.AddEntry(gr6,"14 layers","p")
legend6.Draw()
c5.Draw()
c5.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_14layers.pdf")

# LAYER 15
c6 = TCanvas( 'c6', 'Energy Resolution', 500, 10, 800, 500 )
gr4 = TGraphErrors( n2, x15, y15, 0, e15)
gr4.SetMarkerColor( 2 )
gr4.SetMarkerStyle( 34 )
gr4.SetLineColor(2)
f4 = TF1("f4","sqrt([0]**2/x + [1]**2)",0,400)
f4.SetParameter(0,1)
f4.SetParameter(1,0.5)
f4.SetLineColor( 2 )
result4 = gr4.Fit(f4, "SQR")
gr4.SetTitle( 'Energy Resolution (a and c only)' )
gr4.GetXaxis().SetTitle( 'E [GeV]' )
gr4.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr4.Draw( 'AP' )
legend4 = TLegend(0.3,0.8,0.5,0.9)
legend4.SetHeader("Legend","C")
legend4.AddEntry(gr4,"15 layers","p")
legend4.Draw()
c6.Draw()
c6.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_15layers.pdf")

# LAYER 20
c7 = TCanvas( 'c7', 'Energy Resolution', 500, 10, 800, 500 )
gr5 = TGraphErrors( n, x20, y20, 0, e20)
gr5.SetMarkerColor( 65 )
gr5.SetMarkerStyle( 47 )
gr5.SetLineColor(65)
f5 = TF1("f5","sqrt([0]**2/x + [1]**2)",0,400)
f5.SetParameter(0,1)
f5.SetParameter(1,0.5)
f5.SetLineColor( 65 )
result5 = gr5.Fit(f5, "SQR")
gr5.SetTitle( 'Energy Resolution (a and c only)' )
gr5.GetXaxis().SetTitle( 'E [GeV]' )
gr5.GetYaxis().SetTitle( '\sigma_{E}/E' )
gr5.Draw( 'AP' )
legend5 = TLegend(0.3,0.8,0.5,0.9)
legend5.SetHeader("Legend","C")
legend5.AddEntry(gr5,"20 layers","p")
legend5.Draw()
c7.Draw()
c7.SaveAs("FitGraphs/ACFit/EnergyResolutionFit_AC_25to400GeV_20layers.pdf")