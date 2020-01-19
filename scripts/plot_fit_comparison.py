#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, get_bin_name

import numpy as np

"""@file docstring
Script for plotting comparisons between different types of fits (histogram-based, analytic and hybrid

@author Andres Tiko <andres.tiko@cern.ch>
"""


def make_category_matrix(catRatios, weighted = True):
  b = np.array(catRatios)
  if weighted:
    return (b[:,0], b[:,1])
  else:
    return (b[:,0], b[:,1]/b[:,1])

def make_header():
    latex = ""
    latex += "\documentclass{beamer}\n"
    latex += "\usepackage[utf8]{inputenc}\n"
    latex += "\usepackage{default}\n"
    latex += "\usepackage{graphbox}\n"
    latex += "\\begin{document}\n"
    latex += "\\beamertemplatenavigationsymbolsempty\n"
    #\usepackage[T1]{fontenc}
    return latex

   

def make_footer():
    return """
    \end{document}
    """

def get_results_file(typestring):
    if len(typestring) == 0:
        return "results_cat.txt"
    else:
        return "results_cat_%s.txt" % typestring
        
def get_plot_file(typestring, isPass = True, log = False):
    filename = ""
    if isPass: filename += "pass_"
    else: filename += "fail_"
    if log: filename += "log_"
    filename += "fit_s"
    if len(typestring) > 0:
         filename += "_"+typestring
    filename += ".png"
    return filename

def make_latex(fit1, fit2):
    ratios1 = read_category_ratios(fit1["indir"] + "/" + get_results_file(fit1["type"]))
    ratios2 = read_category_ratios(fit2["indir"] + "/" + get_results_file(fit2["type"]))

    latex = make_header()
    for b in range(21):
        if b == 6: continue
        textcolor = "black"
        if ratios1[b][0] + ratios1[b][1] < ratios2[b][0] - ratios2[b][1] or \
           ratios1[b][0] - ratios1[b][1] > ratios2[b][0] + ratios2[b][1]:
            textcolor = "red"
                
    
        latex += "\\begin{frame}{%s}\n" % get_bin_name(b).replace("_", " ")
        latex += "\\begin{columns}[T,onlytextwidth]\n"
        latex += "\column{0.55\\textwidth}\n"
        latex += "%s:\\\\" % fit1["title"]
        latex += "Pass: \includegraphics[align=c,width=0.8\\textwidth]{%s/bin%d/%s}\\\\ \n" % (fit1["indir"], b, get_plot_file(fit1["type"], True))
        latex += "$ \\textcolor{%s}{%.5f \pm %.5f} $  \\\\ \n" % (textcolor, ratios1[b][0], ratios1[b][1])
        latex += "Fail: \includegraphics[align=c,width=0.8\\textwidth]{%s/bin%d/%s}\\\\ \n" % (fit1["indir"], b, get_plot_file(fit1["type"], False))
        latex += "\column{0.45\\textwidth}\n"
        latex += "%s:\\\\" % fit2["title"]
        latex += "\includegraphics[width=\\textwidth]{%s/bin%d/%s}\\\\ \n" % (fit2["indir"], b, get_plot_file(fit2["type"], True))
        latex += "$ \\textcolor{%s}{%.5f \pm %.5f} $ \n" % (textcolor, ratios2[b][0], ratios2[b][1])
        latex += "\includegraphics[width=\\textwidth]{%s/bin%d/%s}\\\\ \n" % (fit2["indir"], b, get_plot_file(fit2["type"], False))
        latex += "\end{columns}\n"
        latex += "\end{frame}\n"
    latex += make_footer()
    return latex
    


if __name__ == "__main__":
  maindir = "/home/andres/tth/chargeFlip/CMSSW_7_4_7/src/tthAnalysis/ChargeFlipEstimation/bin/"
  #fit1 = {"indir": maindir+"fit_output_data_shiftPeak", "type": "hybrid", "title": "Data hybrid"}
  
  #Input fits to compare  
  fit1 = {"indir": maindir+"fit_output_data_shiftPeak", "type": "shapes", "title": "Analytic fit"}
  fit2 = {"indir": maindir+"fit_output_data_shiftPeak", "type": "", "title": "Histogram-based fit"}
  #Where to write output latex  
  outfile = "latex_output/data_analytic_vs_histos.tex"
  f = open(outfile, "w")
  f.write(make_latex(fit1, fit2))
