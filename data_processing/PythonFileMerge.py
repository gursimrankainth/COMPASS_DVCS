import ROOT
import os
import re

input_dir = "/eos/user/g/gkainth/condorOutput/P06/real/mod/"
chunk_size = 1500

tree_names = ["USR970"] #"USR970" "USR970_gen" "USR971"

files = sorted([
  f for f in os.listdir(input_dir)
  if re.search(r"\.root(\.\d+)?$", f)
])

for tree_name in tree_names:

  print("\nProcessing tree:", tree_name)
  print("Total files found:", len(files))

  if not files:
    print("No ROOT files found, skipping")
    continue

  first_file_path = os.path.join(input_dir, files[0])
  check_chain = ROOT.TChain(tree_name)
  n_added = check_chain.Add(first_file_path)
  if n_added == 0:
    print(f"Tree '{tree_name}' not found in first input file. Stopping.")
    break

  for chunk_idx, start in enumerate(range(0, len(files), chunk_size)):
    chunk_files = files[start:start + chunk_size]
    chain = ROOT.TChain(tree_name)

    for f in chunk_files:
      chain.Add(os.path.join(input_dir, f))

    print(f"Chunk {chunk_idx}: adding {len(chunk_files)} files")
    if tree_name == "USR970_gen":
      output_file = f"merged_gen_chunk_{chunk_idx}.root"
    elif len(tree_names) == 1:
      output_file = f"merged_chunk_{chunk_idx}.root"
    else:
      output_file = f"merged_chunk_{chunk_idx}.root"

    merged_entries = chain.Merge(output_file)
    if merged_entries < 0:
      print(f"Merge failed for {tree_name} chunk {chunk_idx}")
      continue

    print("Wrote:", output_file)
