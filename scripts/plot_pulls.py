#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, bin_names_composite, bin_names_single, \
                                                   mkdir_p, readMisIDRatios
from tthAnalysis.ChargeFlipEstimation.plot_pulls import get_other_component, get_all_containers

import ROOT

"""@file docstring
Script for plotting pulls, superseded by plot_pulls_all.py

@author Andres Tiko <andres.tiko@cern.ch>
"""

def get_bin_nr_composite(cat):
  return bin_names_composite.index(cat)

def get_bin_nr_single(cat):
  return bin_names_single.index(cat)

def make_pull_plot(category, misIDRatios, catRatios, datastring, fitname, fittype):
  pull_plot = ROOT.TH1D(category, category, 6, 0, 6 );
  others_plot = ROOT.TH1D(category+"others", category+"others", 6, 0, 6 );
  bin_names = get_all_containers(category)
  for b in range(1, len(bin_names)+1):
    pull_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    (value, err) = catRatios[get_bin_nr_composite(bin_names[b-1])]
    pull_plot.SetBinContent(b, value)
    pull_plot.SetBinError(b, err)

    other = get_other_component(category, bin_names[b-1])
    (valueO, errO) = misIDRatios[get_bin_nr_single(other)]
    others_plot.SetBinContent(b, valueO)
    others_plot.SetBinError(b, errO)
    #print bin_names[b-1], value, valueO  
  pull_plot.Add(others_plot, -1)
  c = ROOT.TCanvas("Plot", "Plot", 800,600)
  ROOT.gStyle.SetOptStat(0)
  pull_plot.Draw()
  if len(fittype) == 0: fittype = "histograms"
  mydir = "pull_plots/%s/%s/%s/" % (fitname, fittype, datastring)
  mkdir_p(mydir)
  c.SaveAs("%s/%s_pulls.pdf" % (mydir, category))
  c.SaveAs("%s/%s_pulls.png" % (mydir, category))

if __name__ == "__main__":
  for datastring in ["data", "pseudodata"]:
    for fitname in ["summer_May11_noWJ"]:
      for fittype in [""]: #, "shapes", "hybrid"]:
        fittypestring = fittype
        if len(fittype) > 0: fittypestring = "_"+fittype
        file_misId = "fit_output_%s_%s/fit_res%s.root" % (datastring, fitname, fittypestring)
        file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, fitname, fittypestring)
        
        misIDRatios = readMisIdRatios(file_misId)
        catRatios = read_category_ratios(file_cats)

        for bin in bin_names_single:
          make_pull_plot(bin, misIDRatios, catRatios, datastring, fitname, fittype)


