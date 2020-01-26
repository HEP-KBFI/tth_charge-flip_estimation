import ROOT

ROOT.gROOT.SetBatch(True)

import array
import errno    
import os
import argparse

"""@file docstring
Utility functions for charge flip estimation

@author Andres Tiko <andres.tiko@cern.ch>
"""

class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):
  def _split_lines(self, text, width):
    if text.startswith('R|'):
      return text[2:].splitlines()
    return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)

"""Bin names in the form they are in the ntuple"""
BIN_NAMES_COMPOSITE = [
  "BB_LL", "BB_ML", "BB_MM", "BB_HL", "BB_HM", "BB_HH",
  "EE_LL", "EE_ML", "EE_MM", "EE_HL", "EE_HM", "EE_HH",
  "BE_LL", "BE_ML", "EB_ML", "BE_MM", "BE_HL", "EB_HL",
  "BE_HM", "EB_HM", "BE_HH",
]

"""Nicer bin names - historical reason why 2 sets used
   TODO: change the format in ntuples"""
BIN_NAMES_COMPOSITE_NICE = [ '{}{}_{}{}'.format(b[0], b[3], b[1], b[4]) for b in BIN_NAMES_COMPOSITE ]

BIN_NAMES_SINGLE = [ "BL", "BM", "BH", "EL", "EM", "EH" ]

def get_component_cats(nice_name):
  return nice_name.split("_")

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def read_category_ratios(file_cats, exclude_bins = []):
  f = open(file_cats)  
  bins = []
  ratios_dict = {}
  ratios = []
  for line in f.readlines():
    spl = line.split(",")
    bin_nr = int(spl[0])
    #Skip some bins
    bin_name = BIN_NAMES_COMPOSITE_NICE[bin_nr]
    if bin_nr in exclude_bins: continue
    if bin_name in exclude_bins: continue
    if float(spl[2]) == 0: 
      spl[2] = 1e-8
      spl[3] = 1e-8
    ratios.append((float(spl[1]), float(spl[2]), float(spl[3])))
    ratios_dict[bin_name] = (float(spl[1]), float(spl[2]), float(spl[3]) )    
  return (ratios, ratios_dict)

def readMisIDRatios(file_misId):
  f = ROOT.TFile(file_misId)
  misIdHisto = f.Get("chargeMisId")
  ratios = {}
  for etaBin in range(1, misIdHisto.GetNbinsY()+1):
    for ptBin in range(1, misIdHisto.GetNbinsX()+1):
      ratios[get_bin_name_single(etaBin, ptBin)] = (misIdHisto.GetBinContent(ptBin, etaBin), misIdHisto.GetBinError(ptBin, etaBin), misIdHisto.GetBinError(ptBin, etaBin))
  return ratios
  
def get_bin_name_single(bin_nr_eta, bin_nr_pt):
  return BIN_NAMES_SINGLE[(bin_nr_eta - 1) * 3 + (bin_nr_pt - 1)]
  
def get_bin_name(bin_nr):
  return BIN_NAMES_COMPOSITE[bin_nr]

def get_bin_name_nice(bin_nr):
  return BIN_NAMES_COMPOSITE_NICE[bin_nr]
    
def bin_names_to_numbers(bin_names):
  return [ BIN_NAMES_COMPOSITE.index(b) for b in bin_names ]

def get_bin_nr(bin_name_nice):
  return BIN_NAMES_COMPOSITE_NICE.index(bin_name_nice)

def make_title(name):
  title = ""  
  if name == "gen":
    title = "Generator-level"
  elif name == "genRec":
    title = "Generator-level wrt reconstructed p_{T} and #eta"
  elif name == "gen_fit":
    title = "Generator-level, misID rates from solving equations"
  elif name == "gen_fit_exclusions":
    title = "Generator-level, misID rates from solving equations, some categories excluded"
  elif name == "pseudodata":
    title = "Pseudodata"
  elif name == "pseudodata_exclusions":
    title = "Pseudodata, some categories excluded"
  elif name == "Data":
    title = "Pseudodata"
  elif name == "Data_exclusions":
    title = "Pseudodata, some categories excluded"
  return title
  

def fit_results_to_file(rates, uncs, output_filename):
    output_file = ROOT.TFile(output_filename, "recreate")
    output_file.cd()

    binsPt = [ 10, 25, 50, 1000 ]
    binsEta  = [ 0, 1.479, 2.5 ]
    NbinsPt  = len(binsPt) - 1
    NbinsEta = len(binsEta) - 1

    h = ROOT.TH2F(
      "chargeMisId", "chargeMisId;p_{T}(e) [GeV];|#eta(e)|",
      NbinsPt, array.array('d',binsPt),
      NbinsEta, array.array('d',binsEta)
    )
    h.SetStats(0)
    h.SetOption("colz;text")
    for bin_idx in range(len(rates)):
        h.SetBinContent((bin_idx % NbinsPt) + 1, (bin_idx / NbinsPt) + 1, rates[bin_idx])
        h.SetBinError  ((bin_idx % NbinsPt) + 1, (bin_idx / NbinsPt) + 1, uncs[bin_idx])
    output_file.Write()
    output_file.Close()

    print("Wrote file: {}".format(output_filename))
