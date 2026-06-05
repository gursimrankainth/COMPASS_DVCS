import argparse
import glob
import os
import ROOT
import numpy as np
from collections import defaultdict

# ********************************
# Parse arguments
parser = argparse.ArgumentParser(description="Filter tree for generated events.")
parser.add_argument("--period", type=str, required=False, default=None,
                    help="Period string (e.g. 09)")

args = parser.parse_args()
period = args.period if args.period else input("Enter the period (e.g. 09): ").strip()

# ********************************
# Build a TChain from all matching files
input_path = "/eos/user/g/gkainth/BH/mergedFiles/P07_muMinus"
pattern = os.path.join(input_path, "merged_gen_chunk_*.root")
file_list = sorted(glob.glob(pattern))
if not file_list:
  raise FileNotFoundError(f"No files found matching pattern '{pattern}' in {input_path}")

print(f"Adding {len(file_list)} files to chain:")
for f in file_list:
  print("   ", f)

tree_name = "USR970_gen"
chain = ROOT.TChain(tree_name)
for f in file_list:
  chain.Add(f)

n_total = chain.GetEntries()
print(f"Processing {n_total} entries in tree '{tree_name}' from {len(file_list)} files...")

# *******************************
# Prepare output file and trees
output_file = f"gen_P{period}.root"
out_file = ROOT.TFile.Open(output_file, "RECREATE")
out_file.cd()
out_tree = chain.CloneTree(0)

# **********************************
# t' = t - t_min (pg. 140 Johannes' thesis)
def get_tMin(Q2, xbj):
  M_p = 0.93827208816  # GeV

  c1 = np.sqrt(1 + (4 * M_p**2 * xbj**2) / Q2)
  c2 = (2 * M_p**2 * xbj**2) / Q2
  c3 = (2 * M_p**2 * xbj) / Q2

  tMin = - (Q2 / xbj) * ((c2 + 1 - c1) / (c3 + 1 - c1))
  return tMin

# ********************************
# Apply phase space cuts 
for idx in range(n_total):
  chain.GetEntry(idx)
  event = chain # alias for clarity

  if not (2 < event.nu_gen < 270):
   continue

  if not (0.5 < event.Q2_gen < 80):
    continue

  t_prime = event.t_gen - get_tMin(event.Q2_gen, event.xbj_gen)
  if not (-1.2 < t_prime < -0.001):
    continue 

  out_tree.Fill()

# *******************************
# Write to output file
out_tree.Write()
out_file.Close()
print(f"Wrote output to {output_file}")