"""Compare the real-data arrays from Gursimran's and Johannes' code.

G_arrays_P04.pkl stores a dictionary with arrays ordered [t, Q2, nu, phi].
J_arrays_P04.pkl stores eight arrays sequentially, with the real arrays first,
ordered [t, phi, nu, Q2].
"""

import pickle
import argparse
from contextlib import redirect_stdout
import numpy as np


def load_arrays(g_file="G_arrays_P04.pkl", j_file="J_arrays_P04.pkl"):
  with open(g_file, "rb") as f:
    g = pickle.load(f)

  with open(j_file, "rb") as f:
    j_arrays = [pickle.load(f) for _ in range(8)]

  # Convert [t, Q2, nu, phi] to [t, phi, nu, Q2].  The phi conventions
  # run in opposite directions and are offset by three bins:
  # colleague_phi = (3 - our_phi) % 8.
  phi_map = (3 - np.arange(8)) % 8

  g_arrays = {}
  for name in ("real", "BH", "lepPi0", "hepPi0"):
    g_arrays[f"{name}_muPlus"] = np.transpose(
      g[f"{name}_muPlus"], (0, 3, 2, 1)
    )[:, phi_map, :, :]
    g_arrays[f"{name}_muMinus"] = np.transpose(
      g[f"{name}_muMinus"], (0, 3, 2, 1)
    )[:, phi_map, :, :]

  j_names = (
    "real_muPlus", "real_muMinus", "BH_muPlus", "BH_muMinus",
    "lepPi0_muPlus", "lepPi0_muMinus", "hepPi0_muPlus", "hepPi0_muMinus",
  )
  j_arrays = dict(zip(j_names, map(np.asarray, j_arrays)))
  return g_arrays, j_arrays


def compare(name, ours, colleague, rtol=1e-5, atol=1e-8):
  if ours.shape != colleague.shape:
    raise ValueError(f"{name}: shape mismatch: {ours.shape} vs {colleague.shape}")

  diff = ours - colleague
  close = np.isclose(ours, colleague, rtol=rtol, atol=atol)
  nonzero = (ours != 0) | (colleague != 0)
  mismatch = nonzero & ~close

  print(f"\n{name}")
  print("  shape:", ours.shape)
  print("  ours min/max/mean:", np.min(ours), np.max(ours), np.mean(ours))
  print("  J    min/max/mean:", np.min(colleague), np.max(colleague), np.mean(colleague))
  print("  allclose:", np.all(close), f"(rtol={rtol}, atol={atol})")
  print("  nonzero bins (ours/J):", np.count_nonzero(ours), np.count_nonzero(colleague))
  print("  mismatching bins:", np.count_nonzero(mismatch))
  print("  max absolute difference:", np.max(np.abs(diff)))
  print("  mean absolute difference:", np.mean(np.abs(diff)))

  if np.any(mismatch):
    print("  first mismatches [t, phi, nu, Q2]:")
    for index in np.argwhere(mismatch)[:10]:
      idx = tuple(index)
      print(f"    {idx}: ours={ours[idx]:.12g}, J={colleague[idx]:.12g}, diff={diff[idx]:.12g}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--g-file", default="G_arrays_P04.pkl")
    parser.add_argument("--j-file", default="J_arrays_P04.pkl")
    parser.add_argument("--out", default="compare_arrays_P04.out")
    args = parser.parse_args()

    g_arrays, j_arrays = load_arrays(args.g_file, args.j_file)
    with open(args.out, "w") as report, redirect_stdout(report):
        print(f"Comparing {args.g_file} against {args.j_file}")
        for name in ("real", "BH", "lepPi0", "hepPi0"):
            for charge in ("muPlus", "muMinus"):
                key = f"{name}_{charge}"
                compare(key, g_arrays[key], j_arrays[key])
    print(f"Comparison written to '{args.out}'")


if __name__ == "__main__":
  main()
