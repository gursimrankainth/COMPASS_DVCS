import ROOT
import pickle
import os 
import numpy as np
import math
import matplotlib.pyplot as plt
from contextlib import redirect_stdout
from typing import Sequence
from tqdm import tqdm

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import dvcs_constants_2016 as const

from dvcs_constant_scan import * 


# Pickle containing the expensive, period-dependent acceptance calculation.
ACCEPTANCE_CACHE_PATH = Path(__file__).with_name("acceptance_cache.pkl")
SUM_CACHE_PATH = Path(__file__).with_name("sum_cache.pkl")


# *******************************************************************
# *                    *** DATA PREP ***                            *
# *******************************************************************
# **********************************
# Constants 
alpha_em = 0.0072973525693 # electromagnetic fine structure constant 
M_mu = 105.6583755e-3 # GeV/c

# **********************************
# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = (os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root"))

# HEPGEN BH MC Data (Generated data)
hepBH_dir = "/Users/gursimran/cern/2016_data/BH/"
gen_files = (os.path.join(hepBH_dir, "gen_P04_muPlus.root"), os.path.join(hepBH_dir, "gen_P04_muMinus.root"),
             os.path.join(hepBH_dir, "gen_P05_muPlus.root"), os.path.join(hepBH_dir, "gen_P05_muMinus.root"),
             os.path.join(hepBH_dir, "gen_P06_muPlus.root"), os.path.join(hepBH_dir, "gen_P06_muMinus.root"),
             os.path.join(hepBH_dir, "gen_P07_muPlus.root"), os.path.join(hepBH_dir, "gen_P07_muMinus.root"),
             os.path.join(hepBH_dir, "gen_P08_muPlus.root"), os.path.join(hepBH_dir, "gen_P08_muMinus.root"),
             os.path.join(hepBH_dir, "gen_P09_muPlus.root"), os.path.join(hepBH_dir, "gen_P09_muMinus.root"))

# HEPGEN BH MC Data (Reconstructed data)
hepBH_files = (os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P06_muPlus.root"), os.path.join(hepBH_dir, "filtered_P06_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P07_muPlus.root"), os.path.join(hepBH_dir, "filtered_P07_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P08_muPlus.root"), os.path.join(hepBH_dir, "filtered_P08_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P09_muPlus.root"), os.path.join(hepBH_dir, "filtered_P09_muMinus.root"))

# HEPGEN Invisible Pi0 MC Data (Reconstructed data)
hepPi0_dir = "/Users/gursimran/cern/2016_data/HepgenPi0/"
hepPi0_files = (os.path.join(hepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P09_muMinus.root"))

# LEPTO Invisible Pi0 MC Data (Reconstructed data)
lepPi0_dir = "/Users/gursimran/cern/2016_data/LeptoPi0/"
lepPi0_files = (os.path.join(lepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P09_muMinus.root"))

# **********************************
# Total Luminosity 
tot_lum_muPlus = np.sum(const.LUMINOSITY_MUPLUS)
tot_lum_muMinus = np.sum(const.LUMINOSITY_MUMINUS)

# **********************************
# Binning scheme (same as what is used for the final accpetance)
# nu: 4 bins of width 5.5 GeV between 10 and 32 GeV
nu_edges = np.linspace(10, 32, 5)
nu_bins = list(zip(nu_edges[:-1], nu_edges[1:]))
# Q2: 4 bins of width 1 (GeV/c)^2 between 1 and 5
Q2_edges = np.linspace(1, 5, 5)
Q2_bins = list(zip(Q2_edges[:-1], Q2_edges[1:]))
# |t|: 4 bins -> each bin should have roughly the same no. of events 
t_edges = [0.08, 0.136, 0.219, 0.36, 0.64]
t_bins = list(zip(t_edges[:-1], t_edges[1:]))
# phi: 8 bins of width pi/4 rad between -pi and pi
phi_edges = np.linspace(-np.pi, np.pi, 9)
phi_bins = list(zip(phi_edges[:-1], phi_edges[1:]))

# Compute the bin widths
# bin size is uneven in t so it cannot be calculated like this 
delta_nu   = np.diff(nu_edges)       # array of size n_nu_bins
delta_Q2   = np.diff(Q2_edges)       # array of size n_Q2_bins
delta_phi  = np.diff(phi_edges)      # array of size n_phi_bins

# Define the phase space (t, Q2, nu, phi)
phaseSpace = (
  len(t_edges) - 1,
  len(Q2_edges) - 1,
  len(nu_edges) - 1,
  len(phi_edges) - 1
)


# *******************************************************************
# *                 *** HELPER FUNCTIONS ***                        *
# *******************************************************************
# **********************************
# Find all phase space bins (more efficient than handling each axis separately)
# Type annotations used for functions are consistent with Python 3.9+.
def getBin(Q2: float, nu: float, t: float, phi: float) -> tuple[int, int, int, int]:
  Q2_bin = int(math.floor((Q2 - Q2_edges[0]) / (Q2_edges[1] - Q2_edges[0])))
  nu_bin = int(math.floor((nu - nu_edges[0]) / (nu_edges[1] - nu_edges[0])))
  t_bin = np.searchsorted(t_edges, abs(t)) - 1
  phi_bin = int(math.floor((phi - phi_edges[0]) / (phi_edges[1] - phi_edges[0])))

  return Q2_bin, nu_bin, t_bin, phi_bin


# *******************************************************************
# *                    *** ACCEPTANCE ***                           *
# *******************************************************************
# **********************************
# Fill arrays with sum of weights in each bin 
def fill_weights(data_type: str = "gen", period: str = "P04", shape: tuple[int, ...] = phaseSpace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  # Fixed binning structure (Axis order: [t][Q2][nu][phi])
  weights_muPlus     = np.zeros(shape, dtype=np.float64)
  weights_muMinus    = np.zeros(shape, dtype=np.float64)
  weights_muPlus_sq  = np.zeros(shape, dtype=np.float64)
  weights_muMinus_sq = np.zeros(shape, dtype=np.float64)

  n_t, n_Q2, n_nu, n_phi = shape

  # Select tree and files
  if data_type == "gen":
    chain = ROOT.TChain("USR970_gen")
    file_list = gen_files

  elif data_type == "rec":
    chain = ROOT.TChain("USR970_filtered")
    file_list = hepBH_files

  else:
    raise ValueError("data_type must be 'gen' or 'rec'")

  for f in file_list:
    if period in f:
      chain.Add(f)

  # Event loop
  for event in chain:
    charge = event.Q_beam
    weight = event.weight_DVCS

    if data_type == "gen":
      Q2  = event.Q2_gen
      nu  = event.nu_gen
      t   = event.t_gen
      phi = event.phi_gg_gen
    else:
      Q2  = event.Q2_fit
      nu  = event.nu_fit
      t   = event.t_fit
      phi = event.phi_gg_fit

    # Normalize phi to [-pi, pi)
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Bin Q2
    Q2_bin, nu_bin, t_bin, phi_bin = getBin(Q2, nu, t, phi)
    if Q2_bin < 0 or Q2_bin >= n_Q2:
      continue

    # Bin nu
    if nu_bin < 0 or nu_bin >= n_nu:
      continue

    # Bin phi
    if phi_bin < 0 or phi_bin >= n_phi:
      continue

    # Bin |t|
    if t_bin < 0 or t_bin >= n_t:
      continue

    # Fill arrays
    if charge == 1:
      weights_muPlus[t_bin, Q2_bin, nu_bin, phi_bin]     += weight
      weights_muPlus_sq[t_bin, Q2_bin, nu_bin, phi_bin]  += weight * weight

    elif charge == -1:
      weights_muMinus[t_bin, Q2_bin, nu_bin, phi_bin]    += weight
      weights_muMinus_sq[t_bin, Q2_bin, nu_bin, phi_bin] += weight * weight

  return (weights_muPlus, weights_muMinus, weights_muPlus_sq, weights_muMinus_sq)

# **********************************
# Get the acceptance 
def compute_acceptance(rec_muPlus: np.ndarray, rec_muMinus: np.ndarray, gen_muPlus: np.ndarray, gen_muMinus: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
  shape = rec_muPlus.shape

  acc_muPlus  = np.zeros(shape, dtype=np.float64)
  acc_muMinus = np.zeros(shape, dtype=np.float64)

  # muPlus 
  mask_muPlus = gen_muPlus != 0.0
  np.divide(rec_muPlus, gen_muPlus, out=acc_muPlus, where=mask_muPlus)

  # muMinus
  mask_muMinus = gen_muMinus != 0.0
  np.divide(rec_muMinus, gen_muMinus, out=acc_muMinus, where=mask_muMinus)

  return acc_muPlus, acc_muMinus

# **********************************
# Get the acceptance error (returns variance and standard deviation)
def compute_acceptance_error(rec_muPlus_sq: np.ndarray, rec_muMinus_sq: np.ndarray,
                             gen_muPlus_sq: np.ndarray, gen_muMinus_sq: np.ndarray,
                             rec_muPlus: np.ndarray, rec_muMinus: np.ndarray,
                             gen_muPlus: np.ndarray, gen_muMinus: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

  shape = rec_muPlus.shape

  var_muPlus  = np.zeros(shape, dtype=np.float64)
  var_muMinus = np.zeros(shape, dtype=np.float64)

  err_muPlus  = np.zeros(shape, dtype=np.float64)
  err_muMinus = np.zeros(shape, dtype=np.float64)

  for idx in np.ndindex(shape):
    # μ+
    weightRec   = rec_muPlus[idx]
    weightGen   = gen_muPlus[idx]
    weight2Rec  = rec_muPlus_sq[idx]
    weight2Gen  = gen_muPlus_sq[idx]

    if weightGen == 0.0:
      varAcc = 0.0
    else:
      errT1 = (1.0 / weightGen)**2 * weight2Rec
      errT2 = (weightRec / (weightGen**2))**2 * weight2Gen
      varAcc = errT1 + errT2

    var_muPlus[idx] = varAcc
    err_muPlus[idx] = np.sqrt(varAcc)

    # μ-
    weightRec   = rec_muMinus[idx]
    weightGen   = gen_muMinus[idx]
    weight2Rec  = rec_muMinus_sq[idx]
    weight2Gen  = gen_muMinus_sq[idx]

    if weightGen == 0.0:
      varAcc = 0.0
    else:
      errT1 = (1.0 / weightGen)**2 * weight2Rec
      errT2 = (weightRec / (weightGen**2))**2 * weight2Gen
      varAcc = errT1 + errT2

    var_muMinus[idx] = varAcc
    err_muMinus[idx] = np.sqrt(varAcc)

  return var_muPlus, var_muMinus, err_muPlus, err_muMinus


def load_cache(cache_path: Path = ACCEPTANCE_CACHE_PATH) -> dict:
  """Load cached acceptance arrays, returning an empty cache if none exists."""
  if not cache_path.exists():
    return {}

  with cache_path.open("rb") as cache_file:
    cache = pickle.load(cache_file)
  if not isinstance(cache, dict):
    raise ValueError(f"Acceptance cache {cache_path} does not contain a dictionary")
  return cache


def save_acceptance_cache(cache: dict, cache_path: Path = ACCEPTANCE_CACHE_PATH) -> None:
  """Save acceptance arrays for all calculated periods."""
  with cache_path.open("wb") as cache_file:
    pickle.dump(cache, cache_file, protocol=pickle.HIGHEST_PROTOCOL)


def get_cached_acceptance(period: str, cache: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  """Return cached acceptance and variance arrays, checking their shape."""
  required = ("acc_muPlus", "acc_muMinus", "varAcc_muPlus", "varAcc_muMinus")
  if period not in cache or any(key not in cache[period] for key in required):
    raise KeyError(period)

  arrays = tuple(np.asarray(cache[period][key]) for key in required)
  if any(array.shape != phaseSpace for array in arrays):
    raise ValueError(f"Acceptance cache for {period} has an unexpected shape")
  return arrays


def get_cached_sums(period: str, cache: dict) -> tuple[np.ndarray, ...]:
  """Return cached event-sum arrays for a period, checking their shapes."""
  required = (
    "real_muPlus", "real_muMinus", "varD_muPlus", "varD_muMinus",
    "BH_muPlus", "BH_muMinus", "varB_muPlus", "varB_muMinus",
    "lepPi0_muPlus", "lepPi0_muMinus", "varL_muPlus", "varL_muMinus",
    "hepPi0_muPlus", "hepPi0_muMinus", "varH_muPlus", "varH_muMinus",
  )
  if period not in cache or any(key not in cache[period] for key in required):
    raise KeyError(period)

  arrays = tuple(np.asarray(cache[period][key]) for key in required)
  if any(array.shape != phaseSpace for array in arrays):
    raise ValueError(f"Sum cache for {period} has an unexpected shape")
  return arrays

# **********************************
# Check that the acceptance and associated error make sense (test by integrating over phi)
# Plotting funtion is used in test_acceptance
def plot_acceptance_integrated_phi(acceptance_muPlus: np.ndarray, acceptance_err_muPlus: np.ndarray,
                                   acceptance_muMinus: np.ndarray, acceptance_err_muMinus: np.ndarray) -> None:
  nu_bin_centers = 0.5 * (nu_edges[:-1] + nu_edges[1:])
  n_t_bins, n_Q2_bins, n_nu_bins = acceptance_muPlus.shape

  fig, axes = plt.subplots(nrows=n_Q2_bins, ncols=n_t_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_Q2_bins):  # rows (Q2 -> y values)    
    for j in range(n_t_bins):  # columns (t -> x)   
      ax = axes[i, j]

      # acceptance_muPlus and acceptance_muMinus shape: (t, Q2, nu)
      y_muPlus = acceptance_muPlus[j, i]
      yerr_muPlus = acceptance_err_muPlus[j, i]
      y_muMinus = acceptance_muMinus[j, i]
      yerr_muMinus = acceptance_err_muMinus[j, i]

      ax.errorbar(nu_bin_centers, y_muPlus, yerr=yerr_muPlus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='red', markersize=5, ecolor='red', label='μ⁺' if i == 0 and j == 0 else "")
      ax.errorbar(nu_bin_centers, y_muMinus, yerr=yerr_muMinus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='black', markersize=5, ecolor='black', label='μ⁻' if i == 0 and j == 0 else "")

      ax.set_ylim(0, 0.8)
      ax.grid(True, linestyle='--', linewidth=0.5)
      ax.axhline(0, color='gray', linewidth=0.5)

  # Legend
  handles, labels = axes[0, 0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

  # Create a new set of axes for the nu and acceptance scale at the top right corner
  nu_axis = fig.add_axes([0.71, 0.92, 0.185, 0.035])  # [left, bottom, width, height]
  acc_axis = fig.add_axes([0.8, 0.86, 0.11, 0.0825])

  # Setup nu axis
  nu_axis.set_xlim(10, 32)
  nu_axis.set_xticks(nu_edges)
  nu_axis.set_xticklabels([f"{edge:.0f}" for edge in nu_edges], fontsize=12)
  nu_axis.set_yticks([])
  nu_axis.yaxis.set_visible(False)
  nu_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
  nu_axis.xaxis.set_label_position('top')
  nu_axis.set_xlabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
  nu_axis.xaxis.tick_top()
  nu_axis.patch.set_facecolor('none')
  for name, spine in nu_axis.spines.items():
    spine.set_visible(name == 'top')
    if name == 'top':
      spine.set_linewidth(1.0)
      spine.set_color('black')

  # Setup acceptance axis
  acc_axis.set_ylim(0, 0.8)
  acc_axis.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
  acc_axis.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8"], fontsize=12)
  acc_axis.set_xticks([])
  acc_axis.xaxis.set_visible(False)
  acc_axis.set_ylabel("Acceptance", fontsize=14, labelpad=15)
  acc_axis.yaxis.set_label_position('right')
  acc_axis.yaxis.tick_right()
  acc_axis.tick_params(axis='y', direction='in', length=5)
  acc_axis.patch.set_facecolor('none')
  for name, spine in acc_axis.spines.items():
    spine.set_visible(name == 'right')
    if name == 'right':
      spine.set_linewidth(1.0)
      spine.set_color('black')

  # Global axes for t and Q2 
  t_axis = fig.add_axes([0.125, 0.09, 0.775, 0.035])
  Q2_axis = fig.add_axes([0.083, 0.11, 0.034, 0.835])

  # t axis setup
  tick_positions = np.arange(len(t_edges))
  t_axis.set_xlim(0, len(t_edges) - 1)
  t_axis.set_xticks(tick_positions)
  t_axis.set_xticklabels([f"{a}" for a in t_edges], fontsize=12)
  t_axis.set_yticks([])
  t_axis.xaxis.tick_bottom()
  t_axis.set_xlabel(r"$|t|$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  t_axis.patch.set_facecolor('none')
  for name, spine in t_axis.spines.items():
    spine.set_visible(name == 'bottom')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  # Q2 axis setup
  Q2_axis.set_ylim(Q2_edges[-1], Q2_edges[0])
  Q2_axis.set_yticks(Q2_edges)
  Q2_axis.set_yticklabels([f"{int(a)}" for a in Q2_edges], fontsize=12)
  Q2_axis.set_xticks([])
  Q2_axis.yaxis.tick_left()
  Q2_axis.set_ylabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  Q2_axis.patch.set_facecolor('none')
  for name, spine in Q2_axis.spines.items():
    spine.set_visible(name == 'left')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
  plt.savefig("acceptance_integrated_phi.png", dpi=300)

# Integrate over phi and plot the acceptance 
def test_acceptance(rec_muPlus: np.ndarray, rec_muMinus: np.ndarray, rec_muPlus_sq: np.ndarray, rec_muMinus_sq: np.ndarray,
                    gen_muPlus: np.ndarray, gen_muMinus: np.ndarray, gen_muPlus_sq: np.ndarray, gen_muMinus_sq: np.ndarray) -> None:
    
    # Sum over phi (same logic as in acceptance.py)
    rec_muPlus_phiInt  = np.sum(rec_muPlus, axis=3)
    rec_muMinus_phiInt = np.sum(rec_muMinus, axis=3)
    gen_muPlus_phiInt  = np.sum(gen_muPlus, axis=3)
    gen_muMinus_phiInt = np.sum(gen_muMinus, axis=3)

    rec_muPlus_sq_phiInt  = np.sum(rec_muPlus_sq, axis=3)
    rec_muMinus_sq_phiInt = np.sum(rec_muMinus_sq, axis=3)
    gen_muPlus_sq_phiInt  = np.sum(gen_muPlus_sq, axis=3)
    gen_muMinus_sq_phiInt = np.sum(gen_muMinus_sq, axis=3)
    print("Finished integrating over phi for all arrays")

    # Acceptance
    acceptance_muPlus = np.zeros_like(gen_muPlus_phiInt)
    acceptance_muMinus = np.zeros_like(gen_muMinus_phiInt)

    nonzero_muPlus = gen_muPlus_phiInt != 0
    nonzero_muMinus = gen_muMinus_phiInt != 0

    acceptance_muPlus[nonzero_muPlus] = rec_muPlus_phiInt[nonzero_muPlus] / gen_muPlus_phiInt[nonzero_muPlus]
    acceptance_muMinus[nonzero_muMinus] = rec_muMinus_phiInt[nonzero_muMinus] / gen_muMinus_phiInt[nonzero_muMinus]
    print("Finished computing acceptance")

    # Error
    var_muPlus, var_muMinus, err_muPlus, err_muMinus = compute_acceptance_error(
      rec_muPlus_sq_phiInt, rec_muMinus_sq_phiInt,
      gen_muPlus_sq_phiInt, gen_muMinus_sq_phiInt,
      rec_muPlus_phiInt, rec_muMinus_phiInt,
      gen_muPlus_phiInt, gen_muMinus_phiInt
    )
    print("Finished computing acceptance error")

    # Plot
    plot_acceptance_integrated_phi(acceptance_muPlus, err_muPlus, acceptance_muMinus, err_muMinus)
    print("Acceptance plot created successfully!")


# *******************************************************************
# *                    *** CROSS SECTION ***                        *
# *******************************************************************
# **********************************
# Transverse virtual photon flux factor 
def getfluxFac(Q2: float, nu: float, y: float, E_mu: float, xbj: float) -> float:
  c1 = (alpha_em / (2 * np.pi)) * ((1 - xbj) / (Q2 * y * E_mu))
  c2 = y**2 * (1 - (2 * M_mu**2 / Q2))
  c3 = (2 / (1 + Q2 / nu**2)) * (1 - y - Q2 / (4 * E_mu**2))
  virFluxFac = c1 * (c2 + c3)
  return virFluxFac

# **********************************
# Binned sum over the unweighted data for a single period (use for real and LEPTO pi0)
def unweighted_sum(files: Sequence[str], data: str, period: str = "P04") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  muPlus_array = np.zeros(phaseSpace, dtype=np.float64)
  muMinus_array = np.zeros(phaseSpace, dtype=np.float64)
  var_muPlus_array = np.zeros(phaseSpace, dtype=np.float64)
  var_muMinus_array = np.zeros(phaseSpace, dtype=np.float64)

  if data == "real":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "lepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'real' or 'lepPi0'")

  for f in files: 
    if period in f:
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin, nu_bin, t_bin, phi_bin = getBin(Q2, nu, t, phi)
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    if t_bin < 0 or t_bin >= len(t_bins):
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      value = 1/fluxFac
      muPlus_array[i1, i2, i3, i4] += value
      var_muPlus_array[i1, i2, i3, i4] += value**2
    elif charge == -1: 
      value = 1/fluxFac
      muMinus_array[i1, i2, i3, i4] += value
      var_muMinus_array[i1, i2, i3, i4] += value**2

  return muPlus_array, muMinus_array, var_muPlus_array, var_muMinus_array

# **********************************
# Binned sum over the weighted data for a single period (use for HEPGEN BH and HEPGEN pi0)
def weighted_sum(files: Sequence[str], data: str, period: str = "P04") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  muPlus_array = np.zeros(phaseSpace, dtype=np.float64)
  muMinus_array = np.zeros(phaseSpace, dtype=np.float64)
  var_muPlus_array = np.zeros(phaseSpace, dtype=np.float64)
  var_muMinus_array = np.zeros(phaseSpace, dtype=np.float64)

  if data == "hepBH":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "hepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'hepBH' or 'hepPi0'")

  for f in files: 
    if period in f: 
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()
    if data == "hepBH":
      weight = event.weight_PAMBH 
    else: 
      weight = event.weight_all  

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin, nu_bin, t_bin, phi_bin = getBin(Q2, nu, t, phi)
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    if t_bin < 0 or t_bin >= len(t_bins):
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      value = weight/fluxFac
      muPlus_array[i1, i2, i3, i4] += value
      var_muPlus_array[i1, i2, i3, i4] += value**2
    elif charge == -1: 
      value = weight/fluxFac
      muMinus_array[i1, i2, i3, i4] += value
      var_muMinus_array[i1, i2, i3, i4] += value**2

  return muPlus_array, muMinus_array, var_muPlus_array, var_muMinus_array

# **********************************
# Get sum term in 4D cross section (per period)
def get_S(real_sum: np.ndarray, BH_sum: np.ndarray, lepPi0_sum: np.ndarray, hepPi0_sum: np.ndarray,
          period_idx: int, charge: str = "muPlus", systematic: str = "CBH",
          fac: float = 1.0) -> np.ndarray:

  if systematic not in systematicOptions:
    raise ValueError(f"Unknown systematic {systematic!r}; choose from {systematicOptions}")

  cbh_fac = fac if systematic == "CBH" else 1.0
  cpi0_hep_fac = fac if systematic == "CPI0_HEP" else 1.0
  cpi0_lep_fac = fac if systematic == "CPI0_LEP" else 1.0
  r_lepto_fac = fac if systematic == "R_LEPTO" else 1.0

  if charge == "muPlus":
    sum_term = (real_sum - cbh_fac*const.CBH_MUPLUS[period_idx]*BH_sum
                - cpi0_lep_fac*const.CPI0_LEP_MUPLUS[period_idx]*r_lepto_fac*const.R_LEPTO*lepPi0_sum
                - cpi0_hep_fac*const.CPI0_HEP_MUPLUS[period_idx]*(1-r_lepto_fac*const.R_LEPTO)*hepPi0_sum)
  elif charge == "muMinus": 
    sum_term = (real_sum - cbh_fac*const.CBH_MUMINUS[period_idx]*BH_sum
                - cpi0_lep_fac*const.CPI0_LEP_MUMINUS[period_idx]*r_lepto_fac*const.R_LEPTO*lepPi0_sum
                - cpi0_hep_fac*const.CPI0_HEP_MUMINUS[period_idx]*(1-r_lepto_fac*const.R_LEPTO)*hepPi0_sum)
  else:
    raise ValueError('Invalid charge, please use "muPlus" or "muMinus"')
  return sum_term

# **********************************
# Get the mean cross section in (|t|, phi)
def compute_sigma_t_phi(Ncorr_ijkl: np.ndarray) -> np.ndarray:
  n_t, n_Q2, n_nu, n_phi = Ncorr_ijkl.shape
  sigma_t_phi = np.zeros((n_t, n_phi), dtype=np.float64)

  delta_Q2_total = Q2_edges[-1] - Q2_edges[0]
  delta_nu_total = nu_edges[-1] - nu_edges[0]

  for it in range(n_t):
    delta_t = t_edges[it+1] - t_edges[it]

    for iphi in range(n_phi):
      # Sum over Q2 and nu
      S = 0.0
      for iq in range(n_Q2):
        for inu in range(n_nu):
          S += Ncorr_ijkl[it, iq, inu, iphi]

      sigma_t_phi[it, iphi] = (
        S /
        (delta_t *
         delta_Q2_total *
         delta_nu_total *
         delta_phi[iphi])
      )

  return sigma_t_phi

# Integrate over phi (discrete sum over phi bins)
def integrate_over_phi(sigma_t_phi: np.ndarray) -> np.ndarray:
  n_t, n_phi = sigma_t_phi.shape
  sigma_t = np.zeros(n_t, dtype=np.float64)

  for it in range(n_t):
    total = 0.0
    for iphi in range(n_phi):
      total += sigma_t_phi[it, iphi] * delta_phi[iphi]
    sigma_t[it] = total

  return sigma_t


# *******************************************************************
# *                *** CROSS SECTION ERROR ***                      *
# *******************************************************************
# **********************************
# Fill the Var(D_ijkl) or Var(L_ijkl) array per period 
def get_var_S(D_array: np.ndarray, B_array: np.ndarray, L_array: np.ndarray, H_array: np.ndarray,
              period_idx: int, charge: str = "muPlus", systematic: str = "CBH",
              fac: float = 1.0) -> np.ndarray:
  if systematic not in systematicOptions:
    raise ValueError(f"Unknown systematic {systematic!r}; choose from {systematicOptions}")

  cbh_fac = fac if systematic == "CBH" else 1.0
  cpi0_hep_fac = fac if systematic == "CPI0_HEP" else 1.0
  cpi0_lep_fac = fac if systematic == "CPI0_LEP" else 1.0
  r_lepto = fac * const.R_LEPTO if systematic == "R_LEPTO" else const.R_LEPTO

  # the arrays are just the sums, NOT the sqaure root of the sums (no need to sqaure them)
  t1 = D_array
  if charge == "muPlus":
    t2 = (cbh_fac * const.CBH_MUPLUS[period_idx])**2 * B_array
    t3 = (cpi0_lep_fac * const.CPI0_LEP_MUPLUS[period_idx])**2 * r_lepto**2 * L_array
    t4 = (cpi0_hep_fac * const.CPI0_HEP_MUPLUS[period_idx])**2 * (1-r_lepto)**2 * H_array
    sum_term = t1 + t2 + t3 + t4
  elif charge == "muMinus": 
    t2 = (cbh_fac * const.CBH_MUMINUS[period_idx])**2 * B_array
    t3 = (cpi0_lep_fac * const.CPI0_LEP_MUMINUS[period_idx])**2 * r_lepto**2 * L_array
    t4 = (cpi0_hep_fac * const.CPI0_HEP_MUMINUS[period_idx])**2 * (1-r_lepto)**2 * H_array
    sum_term = t1 + t2 + t3 + t4
  else:
    raise ValueError('Invalid charge, please use "muPlus" or "muMinus"')
  return sum_term

# Get variance for 4D cross section per period
def compute_4D_variance(S_ijkl: np.ndarray, A_ijkl: np.ndarray,
                        varS_ijkl: np.ndarray, varA_ijkl: np.ndarray) -> np.ndarray:
    Ncorr = S_ijkl / A_ijkl
    var4D = np.zeros_like(Ncorr)

    # Mask bins with negligible signal
    mask = (S_ijkl != 0) & (A_ijkl != 0)

    # Standard variance propagation formula
    var4D[mask] = (Ncorr[mask]**2) * (
      varS_ijkl[mask] / (S_ijkl[mask]**2) +
      varA_ijkl[mask] / (A_ijkl[mask]**2)
    )

    return var4D

# **********************************
# Get the variance for the mean cross section in |t|, phi 
def get_var_t_phi(var4D_array: np.ndarray) -> np.ndarray:
  n_t, n_Q2, n_nu, n_phi = var4D_array.shape
  var_t_phi = np.zeros((n_t, n_phi), dtype=np.float64)

  delta_Q2_total = Q2_edges[-1] - Q2_edges[0]
  delta_nu_total = nu_edges[-1] - nu_edges[0]

  for it in range(n_t):
    delta_t = t_edges[it+1] - t_edges[it]

    for iphi in range(n_phi):
      # Sum over Q2 and nu
      total_var = 0.0
      for iq in range(n_Q2):
        for inu in range(n_nu):
          total_var += var4D_array[it, iq, inu, iphi]

      norm = delta_t * delta_Q2_total * delta_nu_total * delta_phi[iphi]
      var_t_phi[it, iphi] = total_var / norm**2

  return var_t_phi

# **********************************
# Get the variance for the |t| depedent cross section 
def get_var_t(var_t_phi: np.ndarray) -> np.ndarray:
  n_t, n_phi = var_t_phi.shape
  var_t = np.zeros(n_t)

  for it in range(n_t):
    sum = 0.0

    for iphi in range(n_phi):
      sum += var_t_phi[it, iphi] * delta_phi[iphi]**2

    var_t[it] = sum 

  return var_t

# *******************************************************************
# *                  *** CROSS SECTION WRAPPER ***                  *
# *******************************************************************
# Wrapper function used to detremine the t-dependent cross section and error 
def computeXSec(systematic: str = "CBH", fac: float = 1.0,
                recalcAcc: bool = False,
                recalcSums: bool = False,
                acceptance_cache_path: Path = ACCEPTANCE_CACHE_PATH,
                sum_cache_path: Path = SUM_CACHE_PATH) -> dict:
  if systematic not in systematicOptions:
    raise ValueError(f"Unknown systematic {systematic!r}; choose from {systematicOptions}")
  
  # Dictionary to store the results 
  dvcs_results = {}
  acceptance_cache = {} if recalcAcc else load_cache(acceptance_cache_path)
  cache_updated = False
  sum_cache = {} if recalcSums else load_cache(sum_cache_path)
  sum_cache_updated = False

  # Initialize t-dependent cross-section arrays 
  total_sigma_t_muPlus = np.zeros((4,), dtype=np.float64)
  total_sigma_t_muMinus = np.zeros((4,), dtype=np.float64)
  total_var_t_muPlus = np.zeros((4,), dtype=np.float64)
  total_var_t_muMinus = np.zeros((4,), dtype=np.float64)

  # Loop over the periods
  periods = const.PERIODS
  for idx, period in enumerate(periods):
    print(idx, period)
    print(systematic, fac)
    # ***********************************************
    # *          *** 4D ACCEPTANCE ***              * 
    # ***********************************************
    try:
      acc_muPlus, acc_muMinus, varAcc_muPlus, varAcc_muMinus = get_cached_acceptance(period, acceptance_cache)
      print(f"Using cached acceptance for {period} from {acceptance_cache_path}")
    except (KeyError, ValueError):
      gen_muPlus, gen_muMinus, gen_muPlus_sq, gen_muMinus_sq = fill_weights(data_type="gen", period=period)
      rec_muPlus, rec_muMinus, rec_muPlus_sq, rec_muMinus_sq = fill_weights(data_type="rec", period=period)
      acc_muPlus, acc_muMinus = compute_acceptance(rec_muPlus, rec_muMinus, gen_muPlus, gen_muMinus)
      varAcc_muPlus, varAcc_muMinus, _, _ = compute_acceptance_error(
        rec_muPlus_sq, rec_muMinus_sq, gen_muPlus_sq, gen_muMinus_sq,
        rec_muPlus, rec_muMinus, gen_muPlus, gen_muMinus)
      acceptance_cache[period] = {
        "acc_muPlus": acc_muPlus.copy(), "acc_muMinus": acc_muMinus.copy(),
        "varAcc_muPlus": varAcc_muPlus.copy(), "varAcc_muMinus": varAcc_muMinus.copy(),
      }
      cache_updated = True

    # ***********************************************
    # *         *** 4D CROSS SECTION ***            * 
    # ***********************************************
    try:
      (real_ijkl_muPlus, real_ijkl_muMinus, varD_ijkl_muPlus, varD_ijkl_muMinus,
       BH_ijkl_muPlus, BH_ijkl_muMinus, varB_ijkl_muPlus, varB_ijkl_muMinus,
       lepPi0_ijkl_muPlus, lepPi0_ijkl_muMinus, varL_ijkl_muPlus, varL_ijkl_muMinus,
       hepPi0_ijkl_muPlus, hepPi0_ijkl_muMinus, varH_ijkl_muPlus, varH_ijkl_muMinus) = get_cached_sums(period, sum_cache)
      print(f"Using cached event sums for {period} from {sum_cache_path}")
    except (KeyError, ValueError):
      # These sums depend on the period, but not on the systematic factor.
      real_ijkl_muPlus, real_ijkl_muMinus, varD_ijkl_muPlus, varD_ijkl_muMinus = unweighted_sum(real_files, data="real", period=period)
      BH_ijkl_muPlus, BH_ijkl_muMinus, varB_ijkl_muPlus, varB_ijkl_muMinus = weighted_sum(hepBH_files, data="hepBH", period=period)
      lepPi0_ijkl_muPlus, lepPi0_ijkl_muMinus, varL_ijkl_muPlus, varL_ijkl_muMinus = unweighted_sum(lepPi0_files, data="lepPi0", period=period)
      hepPi0_ijkl_muPlus, hepPi0_ijkl_muMinus, varH_ijkl_muPlus, varH_ijkl_muMinus = weighted_sum(hepPi0_files, data="hepPi0", period=period)
      sum_cache[period] = {
        "real_muPlus": real_ijkl_muPlus.copy(), "real_muMinus": real_ijkl_muMinus.copy(),
        "varD_muPlus": varD_ijkl_muPlus.copy(), "varD_muMinus": varD_ijkl_muMinus.copy(),
        "BH_muPlus": BH_ijkl_muPlus.copy(), "BH_muMinus": BH_ijkl_muMinus.copy(),
        "varB_muPlus": varB_ijkl_muPlus.copy(), "varB_muMinus": varB_ijkl_muMinus.copy(),
        "lepPi0_muPlus": lepPi0_ijkl_muPlus.copy(), "lepPi0_muMinus": lepPi0_ijkl_muMinus.copy(),
        "varL_muPlus": varL_ijkl_muPlus.copy(), "varL_muMinus": varL_ijkl_muMinus.copy(),
        "hepPi0_muPlus": hepPi0_ijkl_muPlus.copy(), "hepPi0_muMinus": hepPi0_ijkl_muMinus.copy(),
        "varH_muPlus": varH_ijkl_muPlus.copy(), "varH_muMinus": varH_ijkl_muMinus.copy(),
      }
      sum_cache_updated = True


    # Acceptance corrected counts
    sum_ijkl_muPlus = get_S(real_ijkl_muPlus, BH_ijkl_muPlus, lepPi0_ijkl_muPlus, 
                                  hepPi0_ijkl_muPlus, period_idx=idx, charge="muPlus", systematic=systematic, fac=fac) 
    Ncorr_ijkl_muPlus = np.zeros_like(sum_ijkl_muPlus)
    np.divide(sum_ijkl_muPlus, acc_muPlus, out=Ncorr_ijkl_muPlus, where=acc_muPlus != 0)

    sum_ijkl_muMinus = get_S(real_ijkl_muMinus, BH_ijkl_muMinus, lepPi0_ijkl_muMinus, 
                                  hepPi0_ijkl_muMinus, period_idx=idx, charge="muMinus", systematic=systematic, fac=fac) 
    Ncorr_ijkl_muMinus = np.zeros_like(sum_ijkl_muMinus)
    np.divide(sum_ijkl_muMinus, acc_muMinus, out=Ncorr_ijkl_muMinus, where=acc_muMinus != 0)

    # ***********************************************
    # *     *** T-DEPDENDENT CROSS SECTION ***      * 
    # ***********************************************
    # mu+ t-dependent cross section 
    sigma_t_phi_muPlus = compute_sigma_t_phi(Ncorr_ijkl_muPlus)
    sigma_t_muPlus = integrate_over_phi(sigma_t_phi_muPlus)
    total_sigma_t_muPlus += sigma_t_muPlus # add per period unnormalized sigma_t to the total for the full 2016 sample 
    sigma_t_muPlus /= const.LUMINOSITY_MUPLUS[idx]
    sigma_t_muPlus *= 1e33 # convert to nb/GeV2
    print(period, "mu+ dsigma/dt (nb/GeV²):", sigma_t_muPlus)

    # mu- t-dependent cross section 
    sigma_t_phi_muMinus = compute_sigma_t_phi(Ncorr_ijkl_muMinus)
    sigma_t_muMinus = integrate_over_phi(sigma_t_phi_muMinus)
    total_sigma_t_muMinus += sigma_t_muMinus # add per period unnormalized sigma_t to the total for the full 2016 sample 
    sigma_t_muMinus /= const.LUMINOSITY_MUMINUS[idx]
    sigma_t_muMinus *= 1e33 # convert to nb/GeV2
    print(period, "mu- dsigma/dt (nb/GeV²):", sigma_t_muMinus)
    
    # ************************************************
    # *       *** 4D CROSS SECTION ERROR ***         * 
    # ************************************************
    varS_ijkl_muPlus = get_var_S(varD_ijkl_muPlus, varB_ijkl_muPlus, varL_ijkl_muPlus, varH_ijkl_muPlus, period_idx=idx, charge="muPlus", systematic=systematic, fac=fac)
    varS_ijkl_muMinus = get_var_S(varD_ijkl_muMinus, varB_ijkl_muMinus, varL_ijkl_muMinus, varH_ijkl_muMinus, period_idx=idx, charge="muMinus", systematic=systematic, fac=fac)

    # ************************************************
    # *   *** T-DEPDENDENT CROSS SECTION ERROR ***   * 
    # ************************************************
    # mu+ t-dependent cross section  error 
    var_ijkl_muPlus = compute_4D_variance(sum_ijkl_muPlus, acc_muPlus, varS_ijkl_muPlus, varAcc_muPlus)
    var_t_phi_muPlus = get_var_t_phi(var_ijkl_muPlus)
    var_t_muPlus = get_var_t(var_t_phi_muPlus)
    total_var_t_muPlus += var_t_muPlus # add per period unnormalized variance_t to the total for the full 2016 sample
    err_t_muPlus = np.sqrt(var_t_muPlus)
    err_t_muPlus /= const.LUMINOSITY_MUPLUS[idx]
    err_t_muPlus *= 1e33 # convert to nb/GeV2
    print(period, "mu+ error (nb/GeV²):", err_t_muPlus)

    # mu- t-dependent cross section error
    var_ijkl_muMinus = compute_4D_variance(sum_ijkl_muMinus, acc_muMinus, varS_ijkl_muMinus, varAcc_muMinus)
    var_t_phi_muMinus = get_var_t_phi(var_ijkl_muMinus)
    var_t_muMinus = get_var_t(var_t_phi_muMinus)
    total_var_t_muMinus += var_t_muMinus # add per period unnormalized variance_t to the total for the full 2016 sample
    err_t_muMinus = np.sqrt(var_t_muMinus)
    err_t_muMinus /= const.LUMINOSITY_MUMINUS[idx]
    err_t_muMinus *= 1e33 # convert to nb/GeV2
    print(period, "mu- error (nb/GeV²):", err_t_muMinus)

    # Save per period results to the dictionary 
    dvcs_results[period] = {
      "sigma_muPlus": sigma_t_muPlus.copy(),
      "sigma_muMinus": sigma_t_muMinus.copy(),
      "err_muPlus": err_t_muPlus.copy(),
      "err_muMinus": err_t_muMinus.copy()
    }

  # ***********************************************
  # *         *** FULL 2016 SAMPLE  ***           * 
  # ***********************************************
  print("Total 2016 Sample")
  # mu +
  total_sigma_t_muPlus /= tot_lum_muPlus
  total_sigma_t_muPlus *= 1e33 # convert to nb/GeV2
  print("total mu+ dsigma/dt (nb/GeV²):", total_sigma_t_muPlus)

  total_err_t_muPlus = np.sqrt(total_var_t_muPlus)
  total_err_t_muPlus /= tot_lum_muPlus
  total_err_t_muPlus *= 1e33 # convert to nb/GeV2
  print("total mu+ error (nb/GeV²):", total_err_t_muPlus)

  # mu-
  total_sigma_t_muMinus /= tot_lum_muMinus
  total_sigma_t_muMinus *= 1e33 # convert to nb/GeV2
  print("total mu- dsigma/dt (nb/GeV²):", total_sigma_t_muMinus)

  total_err_t_muMinus = np.sqrt(total_var_t_muMinus)
  total_err_t_muMinus /= tot_lum_muMinus
  total_err_t_muMinus *= 1e33 # convert to nb/GeV2
  print("total mu- error (nb/GeV²):", total_err_t_muMinus)

  # Save the total 2016 results to the dictionary 
  dvcs_results["total"] = {
    "sigma_muPlus": total_sigma_t_muPlus.copy(),
    "sigma_muMinus": total_sigma_t_muMinus.copy(),
    "err_muPlus": total_err_t_muPlus.copy(),
    "err_muMinus": total_err_t_muMinus.copy()
  }

  if cache_updated:
    save_acceptance_cache(acceptance_cache, acceptance_cache_path)
    print(f"Saved acceptance cache to {acceptance_cache_path}")
  if sum_cache_updated:
    save_acceptance_cache(sum_cache, sum_cache_path)
    print(f"Saved event-sum cache to {sum_cache_path}")

  return dvcs_results


# *******************************************************************
# *                     *** MAIN ***                                *
# *******************************************************************
# **********************************
# Main function - exclude or include functions here 

# Choose which systematic will be studied 
systematicOptions = ("NONE", "CBH", "CPI0_HEP", "CPI0_LEP", "R_LEPTO",)

def main(systematic: str = "CBH", recalcAcc: bool = False, recalcSums: bool = False): 
  if systematic not in systematicOptions:
    raise ValueError(f"Unknown systematic {systematic!r}; choose from {systematicOptions}")
  
  systematic_results = {}
  log_path = f"{systematic}_dvcs_xSection.log"

  with open(log_path, "w") as log_file, redirect_stdout(log_file):
    if systematic == "NONE":
      fac = 1.0
      systematic_results[fac] = computeXSec(systematic=systematic, fac=fac, recalcAcc=recalcAcc, recalcSums=recalcSums)

    elif systematic == "CBH": 
      factors = CBH_FACTORS
      for fac in tqdm(factors, desc=f"{systematic} factors", unit="factor"):
        systematic_results[fac] = computeXSec(systematic=systematic, fac=fac, recalcAcc=recalcAcc, recalcSums=recalcSums)

    else:
      print("Not available yet, script still under developement.")
      return

    results_path = f"{systematic}_dvcs_xSection_results.pkl"
    with open(results_path, "wb") as f:
      pickle.dump(systematic_results, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Results written to '{results_path}'")

  print(f"Log written to '{log_path}'")

if __name__ == "__main__":
  main()
