import argparse
import glob
import os
import ROOT
import numpy as np
from collections import defaultdict
from dataclasses import dataclass

ROOT.ROOT.EnableImplicitMT()

# ********************************
# Parse arguments
parser = argparse.ArgumentParser(description="Filter tree to select DVCS/exclusive photon or pi0 events.")
parser.add_argument(
  "--data", type=str, required=True, choices=["real", "hepgen", "lepto"],
  help="Specify whether the input data is 'real', 'hepgen', or 'lepto'"
)
parser.add_argument("--period", type=str, required=False, default=None,
                    help="Period string (e.g. 09)")

parser.add_argument(
  "--mode", type=str, default="both", choices=["dvcs", "pi0", "both"],
  help="Which tree(s) to produce: 'dvcs', 'pi0', or 'both' (default)"
)

args = parser.parse_args()

# Determine which trees to produce
produce_dvcs = args.mode in ("dvcs", "both")
produce_pi0  = args.mode in ("pi0", "both")

# ********************************
tree_name = "USR970"
period = args.period if args.period else input("Enter the period (e.g. 09): ").strip()

# ********************************
# Build a TChain from all matching files
input_path = "/eos/user/g/gkainth/HepgenPi0/mergedFiles/P09_muMinus/"
pattern = os.path.join(input_path, "merged_chunk_*.root")
file_list = sorted(glob.glob(pattern))
if not file_list:
  raise FileNotFoundError(f"No files found matching pattern '{pattern}' in {input_path}")

print(f"Adding {len(file_list)} files to chain:")
for f in file_list:
  print("   ", f)

chain = ROOT.TChain(tree_name)
for f in file_list:
  chain.Add(f)

n_total = chain.GetEntries()
print(f"Processing {n_total} entries in tree '{tree_name}' from {len(file_list)} files...")

# ********************************
# Determine if real data (set using user input)
is_real_data = args.data == "real"
is_lepto_data = args.data == "lepto"

# ********************************
# Load bad spill list if real data
bad_spills = set()
if is_real_data:
  bad_spill_dir = "/afs/cern.ch/user/g/gkainth/phastPackages/bad_spill"
  pattern = os.path.join(bad_spill_dir, f"P{period}*bad_spill.lst")
  matching_files = glob.glob(pattern)
  if len(matching_files) == 0:
    raise FileNotFoundError(f"No bad spill file found for period P{period} in {bad_spill_dir}")
  elif len(matching_files) > 1:
    print(f"Warning: Multiple bad spill files found for period P{period}. Using the first one: {matching_files[0]}")

  bad_spill_file = matching_files[0]
  print(f"Using bad spill file: {bad_spill_file} ... ")

  with open(bad_spill_file, "r") as f:
    for line in f:
      parts = line.strip().split()
      run, spill = int(parts[0]), int(parts[1])
      flux, LAST, LT, MT, OT, RICH, ECAL, empty = map(int, parts[2:])
      if flux == 1 or LT == 1 or MT == 1 or OT == 1 or ECAL == 1 or empty == 1:
        bad_spills.add((run, spill))

print(f"Loaded {len(bad_spills)} bad spills ...")

# *******************************
# Stats
# triggers, hodoscope, exclusivity (all), fit, Q2, y, t, nu, multiplicity, vis. pi0, excl. photon 
stats_cand = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # counting candidates 
stats_evt = { # counting events 
  "trigger": set(),
  "hodo": set(),
  "excl": set(),
  "fit": set(),
  "Q2": set(),
  "y": set(),
  "t": set(),
  "nu": set(),
  "mult": set(),
  "pi0": set(),
  "dvcs": set()
}

@dataclass(frozen=True)
class Candidate:
  idx: int
  event_key: tuple
  combo_key: tuple

def event_key(event):
  if (period == "09"):
    return (int(event.Run), int(event.Spill), int(event.Evt), str(event.mDST))
  return (int(event.Run), int(event.Spill), int(event.Evt))

def combo_key(candidate):
  return (candidate.inMu_p, candidate.outMu_p, candidate.proton_p, candidate.gamma_p, candidate.gamma_low_p)

# *******************************
# PASS 1: apply all cuts up to multiplicity and record passing candidates
candidates = []
events = defaultdict(list)

for idx in range(n_total):
  chain.GetEntry(idx)
  event = chain  # alias for clarity

  #print("ALL:", event.Evt, event.inMu_TL.P(), event.outMu_TL.P(), event.p_camera_TL.P(), event.gamma_TL.P(), event.gammaLow_TL.P())
  #print("ALL:", event.Evt, event.TiS_flag, event.mDST)

  # Real data checks
  if is_real_data:
    if (int(event.Run), int(event.Spill)) in bad_spills:
      continue

    if not event.TiS_flag:
      continue

  # LEPTO data checks 
  if is_lepto_data:
    if (event.exclLEPTO): # remove exclusive events from the sample 
      continue

  # Trigger + hodoscope
  if not (event.trig_MT or event.trig_OT or event.trig_LT):
    continue
  stats_cand[0] += 1
  stats_evt["trigger"].add(event.Evt)

  if not event.hodoPass:
    continue
  stats_cand[1] += 1
  stats_evt["hodo"].add(event.Evt)

  # Exclusivity variables
  if np.abs(event.delta_pt) > 0.3:
    continue
  delta_phi = event.delta_phi
  if np.abs(delta_phi) > 0.4:
    continue
  if np.abs(event.delta_Z) > 16:
    continue
  if np.abs(event.M2x) > 0.3:
    continue
  stats_cand[2] += 1
  stats_evt["excl"].add(event.Evt)

  # Kinematic fit 
  if not event.fit_conv:
    continue
  chi2_red = event.chi2_fit / event.ndf_fit
  if chi2_red >= 10:
    continue
  stats_cand[3] += 1
  stats_evt["fit"].add(event.Evt)

  if not (1 < event.Q2_fit < 10):
    continue
  stats_cand[4] += 1
  stats_evt["Q2"].add(event.Evt)

  if not (0.05 < event.y_fit < 0.95):
    continue
  stats_cand[5] += 1
  stats_evt["y"].add(event.Evt) 

  if not (0.08 < np.abs(event.t_fit) < 0.64):
    continue
  stats_cand[6] += 1
  stats_evt["t"].add(event.Evt)

  if not (10 < event.nu_fit < 144):
    continue
  stats_cand[7] += 1
  stats_evt["nu"].add(event.Evt)

  #print("PRE:", event.Evt, event.inMu_TL.P(), event.outMu_TL.P(), event.p_camera_TL.P(), event.gamma_TL.P(), event.gammaLow_TL.P())
  # If we reach here, this row passed all pre-multiplicity checks.
  cand = Candidate(
    idx = idx,
    event_key = event_key(event),
    combo_key = (event.inMu_TL.P(), event.outMu_TL.P(), event.p_camera_TL.P(), event.gamma_TL.P(), event.gammaLow_TL.P()),
  )

  candidates.append(cand)
  events[cand.event_key].append(cand)
  
print(f"Pass 1 complete: found {len(candidates)} surviving candidates in {len(events)} events.")


# *******************************
# Prepare output file and trees
output_file = f"filtered_P{period}.root"
out_file = ROOT.TFile.Open(output_file, "RECREATE")

tree_name_dvcs = "USR970_filtered"
tree_name_pi0 = "USR970_pi0"

out_tree_dvcs = chain.CloneTree(0)
out_tree_pi0 = chain.CloneTree(0)

n_kept_dvcs = 0
n_kept_pi0 = 0

# *******************************
# PASS 2: classify surviving candidates by event and fill output trees
def is_exclusive(values):
  return len(set(values)) == 1

saved_events = set()

for event_key, event_candidates in events.items(): # loop over the events 
  inMu_ps = [] 
  outMu_ps = [] 
  proton_ps = [] 
  gamma_ps = [] 

  for candidate in event_candidates: # loop over candidates to check multiplicity 
    # combo_key = candidate.inMu_p, candidate.outMu_p, candidate.proton_p, candidate.gamma_p, candidate.gamma_low_p
    inMu_ps.append(candidate.combo_key[0])
    outMu_ps.append(candidate.combo_key[1])
    proton_ps.append(candidate.combo_key[2])
    gamma_ps.append(candidate.combo_key[3])

  # compute the multiplicity once and apply cut 
  excl_evt = (
    is_exclusive(inMu_ps) and
    is_exclusive(outMu_ps) and
    is_exclusive(proton_ps) and
    is_exclusive(gamma_ps)
  )
   
  if excl_evt == False: 
    continue 

  stats_cand[8] += 1
  stats_evt["mult"].add(event_key[-1])

  #print(event_key, len(event_candidates))

  # check for visible pi0 and save each surviving candidate to the pi0 tree 
  vis_pi0 = False
  for candidate in event_candidates:
    chain.GetEntry(candidate.idx)
    event = chain

    if (event.excl_calo != 2):
      if (event.low_calo == 0 and event.clusterLow_TL.T() > 0.5) or \
          (event.low_calo == 1 and event.clusterLow_TL.T() > 0.63):

        pi0_TL = event.gamma_TL + event.gammaLow_TL
        pi0_M = pi0_TL.M()

        if 0.1061 < pi0_M < 0.1605:
          vis_pi0 = True
          #print("PI0: ", event_key, first_candidate.combo_key)

          if produce_pi0:
            out_tree_pi0.Fill()
            n_kept_pi0 += 1

  # ignore visible pi0 events 
  if vis_pi0:
    continue

  # save DVCS (only first candidate)
  if produce_dvcs and event_key not in saved_events:
    first_candidate = event_candidates[0]
    #print("DVCS:", event_key, first_candidate.combo_key)
    #print(event_key[2])
    chain.GetEntry(first_candidate.idx)
    out_tree_dvcs.Fill()
    n_kept_dvcs += 1
    stats_cand[10] += 1
    stats_evt["dvcs"].add(event.Evt)
    saved_events.add(event_key)
  
print(f"Pass 2 complete: DVCS = {n_kept_dvcs}, (non-exclusive) pi0 = {n_kept_pi0}")


# *******************************
# PASS 3: remove non-exclusive pi0 events
if produce_pi0 and n_kept_pi0 > 0:

  counts = defaultdict(int)

  # prepare reader
  for i in range(out_tree_pi0.GetEntries()):
    out_tree_pi0.GetEntry(i)
    key = (out_tree_pi0.Run, out_tree_pi0.Spill, out_tree_pi0.Evt)
    counts[key] += 1

  out_tree_pi0_excl = out_tree_pi0.CloneTree(0)
  n_kept_pi0_excl = 0

  for i in range(out_tree_pi0.GetEntries()):
    out_tree_pi0.GetEntry(i)
    key = (out_tree_pi0.Run, out_tree_pi0.Spill, out_tree_pi0.Evt)

    if counts[key] != 1:
      continue

    out_tree_pi0_excl.Fill()
    n_kept_pi0_excl += 1
    stats_cand[9] += 1
    stats_evt["pi0"].add(out_tree_pi0.Evt)

  print(f"Pass 3 complete: exclusive pi0 = {n_kept_pi0_excl}")
  out_tree_pi0 = out_tree_pi0_excl

# *******************************
# Write trees
if produce_dvcs:
  out_tree_dvcs.SetName(tree_name_dvcs)
  out_tree_dvcs.Write()

if produce_pi0:
  out_tree_pi0.SetName(tree_name_pi0)
  out_tree_pi0.Write()

out_file.Close()
print(f"Wrote output to {output_file}")

# *******************************
# Print stats
print("\nCutflow:")
print(f"{'Cut':<10} {'Events':>12} {'Candidates':>12}")
for i, cut in enumerate(stats_evt):
  cand = stats_cand[i]
  evt = len(stats_evt[cut])
  print(f"{cut:<10} {evt:>12} {cand:>12}") 