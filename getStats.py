import ROOT
from ROOT import TLorentzVector
import os 
import math
import numpy as np
import matplotlib.pyplot as plt
from dvcs_constants_2016 import *
from scipy.optimize import curve_fit

# **********************************
# Real data 
real_dir = "/eos/user/g/gkainth/real/"
real_files = [os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root")]

# HEPGEN BH MC Data (Generated data)
gen_dir = "/eos/user/g/gkainth/BH/"
gen_files = [os.path.join(gen_dir, "gen_P04_muPlus.root"), os.path.join(gen_dir, "gen_P04_muMinus.root"),
             os.path.join(gen_dir, "gen_P05_muPlus.root"), os.path.join(gen_dir, "gen_P05_muMinus.root"),
             os.path.join(gen_dir, "gen_P06_muPlus.root"), os.path.join(gen_dir, "gen_P06_muMinus.root"),
             os.path.join(gen_dir, "gen_P07_muPlus.root"), os.path.join(gen_dir, "gen_P07_muMinus.root"),
             os.path.join(gen_dir, "gen_P08_muPlus.root"), os.path.join(gen_dir, "gen_P08_muMinus.root"),
             os.path.join(gen_dir, "gen_P09_muPlus.root"), os.path.join(gen_dir, "gen_P09_muMinus.root")]

# HEPGEN BH MC Data (Reconstructed data)
hepBH_dir = "/eos/user/g/gkainth/BH"
hepBH_files = [os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P06_muPlus.root"), os.path.join(hepBH_dir, "filtered_P06_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P07_muPlus.root"), os.path.join(hepBH_dir, "filtered_P07_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P08_muPlus.root"), os.path.join(hepBH_dir, "filtered_P08_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P09_muPlus.root"), os.path.join(hepBH_dir, "filtered_P09_muMinus.root")]

# HEPGEN Invisible Pi0 MC Data (Reconstructed data)
hepPi0_dir = "/eos/user/g/gkainth/HepgenPi0/"
hepPi0_files = [os.path.join(hepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P09_muMinus.root")]

# LEPTO Invisible Pi0 MC Data (Reconstructed data)
lepPi0_dir = "/eos/user/g/gkainth/LeptoPi0/"
lepPi0_files = [os.path.join(lepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P09_muMinus.root")]

# **********************************
# Counting functions 
def count_real_events(nu_bins, period="P04", mode="DVCS"):
  counts_data = [0, 0, 0]
  total_count = 0

  for filename in real_files:
    if period in filename:
      file = ROOT.TFile.Open(filename)
      if mode == "pi0":
        tree = file.Get("USR970_pi0")
      elif mode == "DVCS":
        tree = file.Get("USR970_filtered")

      for event in tree:
        nu = event.nu_fit

        for i, (nu_min, nu_max) in enumerate(nu_bins):
          if nu_min <= nu < nu_max:
            counts_data[i] += 1
            total_count += 1
            break  # avoid double counting

  # Print results
  print(f"\n=== Stats for period {period} ===")
  for i, (nu_min, nu_max) in enumerate(nu_bins):
    print(f"nu in [{nu_min}, {nu_max}): {counts_data[i]} events")

  print(f"Total events: {total_count}")

  return counts_data, total_count

def count_BH_events(nu_bins, period="P04"):
  counts_mc = [0] * len(nu_bins)
  total_count = 0.0

  idx = int(period[1:]) - 4

  for filename in hepBH_files:
    if period not in filename:
      continue

    # Determine charge and corresponding normalization
    if "muPlus" in filename:
      C = CBH_MUPLUS[idx]
    elif "muMinus" in filename:
      C = CBH_MUMINUS[idx]
    else:
      continue  # safety: skip unknown files

    file = ROOT.TFile.Open(filename)
    tree = file.Get("USR970_filtered")

    for event in tree:
      nu = event.nu_fit
      weight = event.weight_PAMBH

      w = weight * C

      for i, (nu_min, nu_max) in enumerate(nu_bins):
        if nu_min <= nu < nu_max:
          counts_mc[i] += w
          total_count += w
          break  # avoid double counting

  # Print results
  print(f"\n=== MC Stats (combined charges) for period {period} ===")
  for i, (nu_min, nu_max) in enumerate(nu_bins):
    print(f"nu in [{nu_min}, {nu_max}): {counts_mc[i]:.6f}")

  print(f"Total normalized MC events: {total_count:.6f}")

  return counts_mc, total_count

def count_HEP_pi0_events(nu_bins, period="P04", mode="DVCS"):
  counts_mc = [0] * len(nu_bins)
  total_count = 0.0

  idx = int(period[1:]) - 4

  for filename in hepPi0_files:
    if period not in filename:
      continue

    # Determine charge and corresponding normalization
    if "muPlus" in filename:
      C = CPI0_HEP_MUPLUS[idx] * (1 - R_LEPTO)
    elif "muMinus" in filename:
      C = CPI0_HEP_MUMINUS[idx] * (1 - R_LEPTO)
    else:
      continue  # safety: skip unknown files

    file = ROOT.TFile.Open(filename)
    if mode == "pi0":
      tree = file.Get("USR970_pi0")
    elif mode == "DVCS":
      tree = file.Get("USR970_filtered")

    for event in tree:
      nu = event.nu_fit
      weight = event.weight_all  

      w = weight * C

      for i, (nu_min, nu_max) in enumerate(nu_bins):
        if nu_min <= nu < nu_max:
          counts_mc[i] += w
          total_count += w
          break  # avoid double counting

  # Print results
  print(f"\n=== MC Stats (combined charges) for period {period} ===")
  for i, (nu_min, nu_max) in enumerate(nu_bins):
    print(f"nu in [{nu_min}, {nu_max}): {counts_mc[i]:.6f}")

  print(f"Total normalized MC events: {total_count:.6f}")

  return counts_mc, total_count

def count_LEP_pi0_events(nu_bins, period="P04", mode="DVCS"):
  counts_mc = [0.0] * len(nu_bins)
  total_count = 0.0

  idx = int(period[1:]) - 4

  for filename in lepPi0_files:
    if period not in filename:
      continue

    # Determine normalization factor once per file
    if "muPlus" in filename:
      C = CPI0_LEP_MUPLUS[idx] * R_LEPTO
    elif "muMinus" in filename:
      C = CPI0_LEP_MUMINUS[idx] * R_LEPTO
    else:
      continue

    file = ROOT.TFile.Open(filename)
    if mode == "pi0":
      tree = file.Get("USR970_pi0")
    elif mode == "DVCS":
      tree = file.Get("USR970_filtered")

    for event in tree:
      nu = event.nu_fit

      # Find bin and accumulate
      for i, (nu_min, nu_max) in enumerate(nu_bins):
        if nu_min <= nu < nu_max:
          counts_mc[i] += C
          total_count += C
          break

  # Output
  print(f"\n=== LEP π0 MC Stats (combined charges) for period {period} ===")
  for i, (nu_min, nu_max) in enumerate(nu_bins):
    print(f"nu in [{nu_min}, {nu_max}): {counts_mc[i]:.6f}")

  print(f"Total normalized MC events: {total_count:.6f}")

  return counts_mc, total_count

# **********************************
# Main loop 
def main():
  nu_bins = [(80,144), (32,80), (10,32)]
  for period in PERIODS:
    # mode can either be DVCS or pi0
    #count_LEP_pi0_events(nu_bins, period=period, mode="DVCS")
    count_HEP_pi0_events(nu_bins, period=period, mode="DVCS")
    #count_BH_events(nu_bins, period=period)
    #count_real_events(nu_bins, period=period, mode="DVCS")

if __name__ == "__main__":
  main()