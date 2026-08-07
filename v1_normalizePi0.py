import os
import glob
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import datetime
from dvcs_constants_2016 import *

# **********************************
# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = (os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root"))

# HEPGEN Pi0 MC Data (Reconstructed data)
hep_dir = "/Users/gursimran/cern/2016_data/HepgenPi0/"
hep_files = (os.path.join(hep_dir, "filtered_P04_muPlus.root"), os.path.join(hep_dir, "filtered_P04_muMinus.root"),
             os.path.join(hep_dir, "filtered_P05_muPlus.root"), os.path.join(hep_dir, "filtered_P05_muMinus.root"),
             os.path.join(hep_dir, "filtered_P06_muPlus.root"), os.path.join(hep_dir, "filtered_P06_muMinus.root"),
             os.path.join(hep_dir, "filtered_P07_muPlus.root"), os.path.join(hep_dir, "filtered_P07_muMinus.root"),
             os.path.join(hep_dir, "filtered_P08_muPlus.root"), os.path.join(hep_dir, "filtered_P08_muMinus.root"),
             os.path.join(hep_dir, "filtered_P09_muPlus.root"), os.path.join(hep_dir, "filtered_P09_muMinus.root"))

# LEPTO Visible Pi0 MC Data (Reconstructed data)
lepVis_dir = "/Users/gursimran/cern/2016_data/LeptoPi0/"
lepVis_files = [os.path.join(lepVis_dir, "filtered_P04_muPlus.root"), os.path.join(lepVis_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P05_muPlus.root"), os.path.join(lepVis_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P06_muPlus.root"), os.path.join(lepVis_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P07_muPlus.root"), os.path.join(lepVis_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P08_muPlus.root"), os.path.join(lepVis_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P09_muPlus.root"), os.path.join(lepVis_dir, "filtered_P09_muMinus.root")]


# **********************************
# ***     c_HEP and c_LEP        ***
# **********************************
# **************
# Get c_HEP 
def hepgenNorm():
  real_mup = []
  hep_mup = []
  real_mum = []
  hep_mum = []
  
  # loop over the real files and count the number of visible pi0's per period 
  for fpath in real_files:
    pos_count = 0
    neg_count = 0

    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      print(f"Failed to open {fpath}")
      continue

    tree = f.Get("USR970_pi0")  

    if not tree:
      print(f"No tree found in {fpath}")
      f.Close()
      continue

    # Example processing
    n_entries = tree.GetEntries()
    for i in range(n_entries):
      tree.GetEntry(i)

      if tree.Q_beam == 1:
        pos_count += 1
      else:
        neg_count += 1

    f.Close()

    real_mup.append(pos_count)
    real_mum.append(neg_count)

  # loop over the hepgen files and count the number of visible pi0's per period
  for fpath in hep_files:
    count = 0

    if "muPlus" in fpath:
      target = hep_mup
    elif "muMinus" in fpath:
      target = hep_mum
    else:
      print("Unknown charge:", fpath)
      continue

    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      continue

    tree = f.Get("USR970_pi0")
    if not tree:
      f.Close()
      continue

    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      count += tree.weight_all

    f.Close()

    target.append(count)

  c_hep_mup = np.array(real_mup) / np.array(hep_mup)
  c_hep_mum = np.array(real_mum) / np.array(hep_mum)
  print(c_hep_mup)
  print(c_hep_mum)
  return c_hep_mup, c_hep_mum

# **************
# Get c_LEP 
def leptoNorm():
  real_mup = []
  lep_mup = []
  real_mum = []
  lep_mum = []
  
  # loop over the real files and count the number of visible pi0's per period 
  for fpath in real_files:
    pos_count = 0
    neg_count = 0

    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      print(f"Failed to open {fpath}")
      continue

    tree = f.Get("USR970_pi0")  
    if not tree:
      print(f"No tree found in {fpath}")
      f.Close()
      continue

    # Example processing
    n_entries = tree.GetEntries()
    for i in range(n_entries):
      tree.GetEntry(i)

      if tree.Q_beam == 1:
        pos_count += 1
      else:
        neg_count += 1

    f.Close()

    real_mup.append(pos_count)
    real_mum.append(neg_count)

  # loop over the lepto files and count the number of visible pi0's per period
  for fpath in lepVis_files:
    count = 0

    if "muPlus" in fpath:
      target = lep_mup
    elif "muMinus" in fpath:
      target = lep_mum
    else:
      print("Unknown charge:", fpath)
      continue

    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      continue

    tree = f.Get("USR970_pi0")
    if not tree:
      f.Close()
      continue

    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      count += 1

    f.Close()

    target.append(count)

  c_lep_mup = np.array(real_mup) / np.array(lep_mup)
  c_lep_mum = np.array(real_mum) / np.array(lep_mum)
  print(c_lep_mum)
  print(c_lep_mup)
  return c_lep_mup, c_lep_mum

# **********************************
#  Write calculated values to Python file for storage
def write_constants_pi0(const_hep_mup, const_lep_mup, const_hep_mum, const_lep_mum, r_lepto, year=2016):
  constants_file = f"/Users/gursimran/cern/COMPASS_DVCS/dvcs_constants_{year}.py" 
  if not os.path.exists(constants_file):
    print(f"No constants file found for {year}: {constants_file}")
    return

  source_file = os.path.basename(__file__)
  today = datetime.date.today().isoformat()  # "YYYY-MM-DD"

  with open(constants_file, "r") as f:
    lines = f.readlines()
  
  # Safety checks
  assert "CBH_MUPLUS" in lines[8]
  assert "CBH_MUMINUS" in lines[9]

  lines[11] = f'SOURCE_FILE = "{source_file} [{today}]"\n'
  lines[12] = f"CPI0_HEP_MUPLUS = {const_hep_mup.tolist()}\n"
  lines[13] = f"CPI0_HEP_MUMINUS = {const_hep_mum.tolist()}\n"
  lines[14] = f"CPI0_LEP_MUPLUS = {const_lep_mup.tolist()}\n"
  lines[15] = f"CPI0_LEP_MUMINUS = {const_lep_mum.tolist()}\n"
  lines[16] = f"R_LEPTO = {r_lepto}\n"

  with open(constants_file, "w") as f:
    f.writelines(lines)

  print("Values written to file:", constants_file)


# **********************************
# ***     Phase space check      ***
# **********************************
def restrictPhaseSpace(q2_min=1, q2_max=10, nu_min=10, nu_max=144):
  """Calculate the pi0 normalization constants in a restricted phase space.

  This is a stand-alone side study: it does not update the constants file and
  is not called by ``main``.  The limits are exclusive, matching the standard
  phase-space definition.

  Returns
  -------
  tuple
    ``(c_hep_mup, c_hep_mum, c_lep_mup, c_lep_mum)``; each item is a
    NumPy array containing one normalization constant per data-taking period.
  """
  if q2_min >= q2_max:
    raise ValueError("q2_min must be smaller than q2_max")
  if nu_min >= nu_max:
    raise ValueError("nu_min must be smaller than nu_max")

  def in_phase_space(event):
    return (q2_min < event.Q2_fit < q2_max and
            nu_min < event.nu_fit < nu_max)

  def get_tree(fpath):
    root_file = ROOT.TFile.Open(fpath)
    if not root_file or root_file.IsZombie():
      raise OSError(f"Failed to open {fpath}")

    tree = root_file.Get("USR970_pi0")
    if not tree:
      root_file.Close()
      raise RuntimeError(f"No USR970_pi0 tree found in {fpath}")

    return root_file, tree

  real_mup = []
  real_mum = []
  for fpath in real_files:
    pos_count = 0
    neg_count = 0
    root_file, tree = get_tree(fpath)

    for event in tree:
      if not in_phase_space(event):
        continue
      if event.Q_beam == 1:
        pos_count += 1
      else:
        neg_count += 1

    root_file.Close()
    real_mup.append(pos_count)
    real_mum.append(neg_count)

  def count_mc(files, weighted):
    counts_mup = []
    counts_mum = []

    for fpath in files:
      count = 0.0
      root_file, tree = get_tree(fpath)

      for event in tree:
        if in_phase_space(event):
          count += event.weight_all if weighted else 1

      root_file.Close()
      if "muPlus" in fpath:
        counts_mup.append(count)
      elif "muMinus" in fpath:
        counts_mum.append(count)
      else:
        raise ValueError(f"Cannot determine beam charge from {fpath}")

    return np.asarray(counts_mup), np.asarray(counts_mum)

  hep_mup, hep_mum = count_mc(hep_files, weighted=True)
  lep_mup, lep_mum = count_mc(lepVis_files, weighted=False)
  real_mup = np.asarray(real_mup)
  real_mum = np.asarray(real_mum)

  with np.errstate(divide="ignore", invalid="ignore"):
    c_hep_mup = real_mup / hep_mup
    c_hep_mum = real_mum / hep_mum
    c_lep_mup = real_mup / lep_mup
    c_lep_mum = real_mum / lep_mum

  print(f"Phase space: {q2_min} < Q2 < {q2_max}, "
        f"{nu_min} < nu < {nu_max}")
  print("c_HEP mu+:", c_hep_mup)
  print("c_HEP mu-:", c_hep_mum)
  print("c_LEP mu+:", c_lep_mup)
  print("c_LEP mu-:", c_lep_mum)

  return c_hep_mup, c_hep_mum, c_lep_mup, c_lep_mum


# **********************************
# ***           Main             ***
# **********************************
# **********************************
# Run the functions and write the constants to the output file
def main():
  c_hep_mup, c_hep_mum = hepgenNorm()
  c_lep_mup, c_lep_mum = leptoNorm()
  r_lepto = 0.4

  #c_hep_mup, c_hep_mum, c_lep_mup, c_lep_mum = restrictPhaseSpace(q2_min=1, q2_max=5, nu_min=10, nu_max=144) 

  write_constants_pi0(const_hep_mup=c_hep_mup, const_hep_mum=c_hep_mum,
                      const_lep_mup=c_lep_mup, const_lep_mum=c_lep_mum,
                      r_lepto=r_lepto)



if __name__ == "__main__":
  main()
