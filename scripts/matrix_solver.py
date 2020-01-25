#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, bin_names_to_numbers, fit_results_to_file
from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_M, make_category_matrix, calculate_rates, \
                                                           make_coefficient_matrix, calculate, calculate_solution

import numpy as np
import math

"""@file docstring
Script for solving overconstrained equation system for charge flip estimation

#See: <https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>
for notation and background

@author Andres Tiko <andres.tiko@cern.ch>
"""

#obsolete
def gen_level_probs():
  return [(0.000035157, 0.000033635), (0.000221524, 0.000014054), (0.000267707, 0.000035297), (0.000829932, 0.000155612), (0.001864478, 0.000074103), (0.002657734, 0.000203647)]

#Make category ratios from gen ratios
#TODO: could be updated to use the function in plot_pulls_gen (akeCatRatiosFrom6)
def make_cat_ratios_from_gen():
  gen_probs = gen_level_probs()
  indices = [[0,0], [0,1], [1,1], [2,0], [2,1], [2,2], [3,3], [4,3], [4,4], [5,3], [5,4], [5,5],
      [0,3], [1,3], [4,0], [1,4], [2,3], [5,0], [2,4], [5,1], [2,5]]
  cat_ratios = []
  for r in range(len(indices)):
    val = gen_probs[indices[r][0]][0]+gen_probs[indices[r][1]][0]
    err = math.sqrt(gen_probs[indices[r][0]][1]**2+gen_probs[indices[r][1]][1]**2)
    if err == 0: err = 1e-8
    cat_ratios.append((val, err))
  return cat_ratios

#Seems to do the same as calculate...cannot remember why two different functions used
def solve_matrix(catRatios):
  A = make_coefficient_matrix()
  (b, W) = make_category_matrix(catRatios)
  W = np.sqrt(np.diag(1/W))
  Aw = np.dot(W,A)
  Bw = np.dot(b,W)
  x_lstsq = np.linalg.lstsq(Aw,Bw)[0] # computing the numpy solution

  Q,R = np.linalg.qr(Aw) # qr decomposition of A
  Qb = np.dot(Q.T,Bw) # computing Q^T*b (project b onto the range of A)
  x_qr = np.linalg.solve(R,Qb) # solving R*x = Q^T*b

  print 'qr solution'
  print x_qr
  print 'lstqs solution'
  print x_lstsq
  return x_qr
  
def print_solution(x):
  for x, value in np.ndenumerate(x):
    print x,value*100

def make_gen_check():
  x = solve_matrix(make_cat_ratios_from_gen())
  print_solution(x)

def print_solution_with_uncertainties(x, uncs):
  for i in range(len(x)):
    print x[i]*100, "+-", uncs[i]*100
    
def print_latex_header():
   print """\\begingroup\setlength{\\fboxsep}{0pt}
    \colorbox{cyan}{%
    \\begin{tabular}{llccc}
	    \hline
	    ID & & $10\leq\pt<25\GeV$ & $25\leq\pt<50\GeV$ & $50\GeV\leq\pt$ \\\\
	    \hline	    
    \end{tabular}%
}\endgroup"""

if __name__ == "__main__":
  #Now run from plot_pulls_all.py
  #Only use this for testing
  exclude_bin_names = []
  exclude_bins = bin_names_to_numbers(exclude_bin_names)
  for datastring in ["pseudodata", "data"]:
    for fitname in ["summer_June6"]:
      for fittype in [""]: #["hybrid", "shapes"] can also be used when fits available
        fittypestring = fittype
        if len(fittype) > 0: fittypestring = "_"+fittype
        file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, fitname, fittypestring)
        categoryRatios = read_category_ratios(file_cats, exclude_bins)
        calculate_solution(categoryRatios, exclude_bins, fitname, fittypestring, datastring)
  print_latex_header()     
