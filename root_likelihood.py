import ROOT
import numpy as np
import array as arr
import pickle
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from dvcs_constants_2016 import * 

# **********************************
# Slope fit function - binned likelihood fit using ROOT 
# ********
# Models  
def log_likelihood(B, N, t_edges):
  t_min = t_edges[0]
  t_max = t_edges[-1]
  denom = np.exp(-B*t_min) - np.exp(-B*t_max)

  mu = np.zeros(len(t_edges) - 1)
  for i in range(len(t_edges) - 1):
    numerator = np.exp(-B * t_edges[i]) - np.exp(-B * t_edges[i+1])
    mu[i] = numerator / denom

  mu = np.clip(mu, 1e-12, None)

  return -np.dot(N, np.log(mu))

def log_likelihood_var(B, sigma_t_error, t_edges):
  t_min = t_edges[0]
  t_max = t_edges[-1]
  denom = np.exp(-B*t_min) - np.exp(-B*t_max)
  norm_factor = sigma_t_error.sum() / denom

  mu = np.zeros(len(t_edges) - 1)
  for i in range(len(t_edges) - 1):
    numerator = np.exp(-B * t_edges[i]) - np.exp(-B * t_edges[i+1])
    mu[i] = norm_factor * numerator

  mu = np.clip(mu, 1e-12, None)

  return -np.dot(sigma_t_error, np.log(mu))

# ********
# ROOT wrappers
class LikelihoodSlope(ROOT.Math.IMultiGenFunction):
  def __init__(self, N, t_edges):
    super().__init__()
    self.N = N
    self.t_edges = t_edges

  def NDim(self):
    return 1

  def Clone(self):
    return LikelihoodSlope(self.N, self.t_edges)

  def DoEval(self, x):
    return log_likelihood(x[0], self.N, self.t_edges)

class LikelihoodVar(ROOT.Math.IMultiGenFunction):
  def __init__(self, sigma_t_error, t_edges):
    super().__init__()
    self.sigma_t_error = sigma_t_error
    self.t_edges = t_edges

  def NDim(self):
    return 1

  def Clone(self):
    return LikelihoodVar(self.sigma_t_error, self.t_edges)

  def DoEval(self, x):
    return log_likelihood_var(x[0], self.sigma_t_error, self.t_edges)

# ********
# Fit functions 
def fit_slope_root(t_edges, dsigma_dt):
  N = np.zeros(len(dsigma_dt))

  for i in range(len(dsigma_dt)):
    delta_t = t_edges[i+1] - t_edges[i]
    N[i] = dsigma_dt[i] * delta_t

  func = LikelihoodSlope(N, t_edges)

  fitter = ROOT.Fit.Fitter()
  params = arr.array('d', [1.0])

  fitter.FitFCN(func, params)

  minimizer = fitter.GetMinimizer()
  #minimizer.SetErrorDef(0.5)
  minimizer.SetErrorDef(0.5)

  minimizer.SetVariable(0, "B", 0.2, 0.001)
  minimizer.Minimize()
  minimizer.Hesse()

  B_fit = minimizer.X()[0]
  cov   = minimizer.CovMatrix(0, 0)

  return B_fit, cov

def fit_variance_root(t_edges, dsigma_dt_err):
  sigma_t_error = np.zeros(len(dsigma_dt_err))

  for i in range(len(dsigma_dt_err)):
    delta_t = t_edges[i+1] - t_edges[i]
    sigma_t_error[i] = dsigma_dt_err[i] * delta_t**2

  func = LikelihoodVar(sigma_t_error, t_edges)

  fitter = ROOT.Fit.Fitter()
  params = arr.array('d', [1.0])

  fitter.FitFCN(func, params)

  minimizer = fitter.GetMinimizer()
  minimizer.SetErrorDef(0.5)

  minimizer.SetVariable(0, "B", 0.2, 0.001)
  minimizer.Minimize()
  minimizer.Hesse()

  B_fit = minimizer.X()[0]
  cov   = minimizer.CovMatrix(0, 0)

  return B_fit, cov

# ********
# Final result 
def dvcs_slope_likelihood_root(t_edges, dsigma_dt, dsigma_dt_err):
  B_fit, cov_main = fit_slope_root(t_edges, dsigma_dt)
  B_var_fit, cov_var = fit_variance_root(t_edges, dsigma_dt_err)

  B_err_combined = np.sqrt(cov_main**2 / cov_var)

  print("B fit =", B_fit)
  print("B combined err =", B_err_combined)

  return B_fit, B_var_fit, B_err_combined