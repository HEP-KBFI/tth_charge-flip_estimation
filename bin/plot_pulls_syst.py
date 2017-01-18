import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT
from utils import read_category_ratios, bin_names_composite as bin_names, bin_names_single, mkdir_p

def make_pull_plots(fittype, syst):
  syst_name = syst.split("_")[2]
  pull_plot = TH1D(syst_name, syst_name, 21, 0, 21 )
  for bin in [0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,20]:
    cat_name = bin_names[bin]
    filename = "%s/bin%d/mlfit.root" % (fittype, bin)
    f = TFile(filename)
    tree = f.Get("tree_fit_sb")
    print tree
    print "%s" % syst
    tree.Draw("%s>>hist"  % syst)
    hist = ROOT.gDirectory.Get("hist")
    print hist
    
    value = hist.GetMean()
    print value
    #print value.GetTree()
    #print value.GetTree().GetHistogram()
    pull_plot.GetXaxis().SetBinLabel(bin + 1, bin_names[bin+1])
    pull_plot.SetBinContent(bin + 1, value)
    f.Close()

  c = TCanvas("Plot", "Plot", 800,600)
  ROOT.gStyle.SetOptStat(0)
  pull_plot.Draw()
  c.SaveAs("plot_pulls/syst_pulls_%s.pdf" % syst_name)
  c.SaveAs("plot_pulls/syst_pulls_%s.png" % syst_name)  
   


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  fittype = "fit_output_data_tightCharge"
  for syst in ["CMS_ttHl_electronER", "CMS_ttHl_electronESBarrel", "CMS_ttHl_electronESEndcap"]:
    make_pull_plots(fittype, syst)
    
