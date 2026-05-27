import ROOT
import os
import glob

# -----------------------------
# USER CONFIG
# -----------------------------

directory = "/eos/user/g/gkainth/condorOutput/P07/real/"
output_directory = "/eos/user/g/gkainth/condorOutput/P07/real/mod"

files = [
  "tree-274509-8-8.root",
  # "tree-272342-81-8.root.001",
]

# -----------------------------
# PROCESS FUNCTION
# -----------------------------

def process_file(input_file):
  base_name = os.path.basename(input_file)

  if not base_name.startswith("tree-"):
    print("Skipping (bad format):", base_name)
    return

  core = base_name[len("tree-"):]

  root_pos = core.find(".root")
  if root_pos == -1:
    print("Skipping (no .root found):", base_name)
    return

  output_file = (
    os.path.join(output_directory, "tree-") +
    core[:root_pos] +
    "_mod" +
    core[root_pos:]
  )

  mdst_value = "mDST-" + core
  tree_name = "USR970"

  fin = ROOT.TFile.Open(input_file)
  if not fin:
    print("Failed to open:", input_file)
    return

  tree = fin.Get(tree_name)
  if not tree:
    print("Tree not found in:", input_file)
    fin.Close()
    return

  fout = ROOT.TFile.Open(output_file, "RECREATE")

  new_tree = tree.CloneTree(0)

  mdst_var = ROOT.std.string()
  new_tree.Branch("mDST", mdst_var)

  for event in tree:
    mdst_var.clear()
    mdst_var.append(mdst_value)
    new_tree.Fill()

  new_tree.Write()
  fout.Close()
  fin.Close()

  print("Processed:", base_name, "->", os.path.basename(output_file))


# -----------------------------
# FILE SELECTION LOGIC
# -----------------------------

if len(files) == 0:
  files_to_process = sorted(
    glob.glob(os.path.join(directory, "tree-*.root*"))
  )
else:
  files_to_process = [
    f if os.path.isabs(f) else os.path.join(directory, f)
    for f in files
  ]

os.makedirs(output_directory, exist_ok=True)


# -----------------------------
# RUN
# -----------------------------

for f in files_to_process:
  process_file(f)
