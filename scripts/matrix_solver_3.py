#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_M, make_category_matrix, calculate_rates

import numpy as np
import math

def make_coefficient_matrix(exclude_bins=[]):
 coeffs = [[ 2,  0],
 [ 1,  1],
 [ 0,  2]]
 used_coeffs = []
 for i in range(len(coeffs)):
   if not i in exclude_bins:
     used_coeffs.append(coeffs[i]) 
 return np.array(used_coeffs)

def make_cat_ratios_from_gen():
  cat_ratios = [(0.0049654093207578565, 0.0001462259953501589),
        (0.0023022851031497876, 0.000041019420446535004), 
        (0.0004516881846517407, 0.000012535791157441956)]
  return cat_ratios

def calculate_uncertainties(M, deltab):
  uncs = [] 
  for i in range(2):    
    uncs.append(0)
    for k in range(len(M[0])):
      uncs[i] += (M[i, k] * deltab[k]) ** 2
    uncs[i] = math.sqrt(uncs[i])
  return np.array(uncs)

def print_solution_latex(x, uncs, datastring):
  latex = """
	    \multirow{2}{*}{%s}   & $0\leq\eta<1.479$    & %.4f $\pm$ %.8f   \\\\
	                            & $1.479\leq\eta<2.5$  & %.4f $\pm$ %.8f   \\\\
	    \hline""" % (datastring, x[0]*100, uncs[0]*100, x[1]*100, uncs[1]*100)
  print latex

def calculate(catRatios, exclude_bins = [], weighted = True):
  A = make_coefficient_matrix(exclude_bins)
  (b, W) = make_category_matrix(catRatios, weighted)
  w = np.diag(1/W**2)
  Aw = np.dot(w, A)  
  M = calculate_M(Aw)
  bw = np.dot(w, b)
  rates = calculate_rates(M, bw)
  (b, err) = make_category_matrix(catRatios)
  uncs = calculate_uncertainties(np.dot(M,w), W)
  return (rates.tolist(), uncs.tolist())


def calculate_solution(categoryRatios, exclude_bins, fitname, fittypestring, datastring):
  (rates, uncs) = calculate(categoryRatios, exclude_bins)
  print_solution_latex(rates, uncs, datastring)   

if __name__ == "__main__":
  exclude_bins = []
  categoryRatios = make_cat_ratios_from_gen()
  calculate_solution(categoryRatios, [], "asd", "asd", "gencheck")
