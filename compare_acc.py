import numpy as np
import pickle 

def load_colleague_acc(filename):
  phi_max = 0
  t_max = 0
  nu_max = 0
  Q2_max = 0

  entries = []

  with open(filename, 'r') as f:
    for line in f:
      if line.startswith('phi') or line.strip() == '':
        continue

      parts = line.split()

      phi_bin = int(parts[0])
      t_bin   = int(parts[1])
      nu_bin  = int(parts[2])
      Q2_bin  = int(parts[3])
      value   = float(parts[4])

      entries.append((Q2_bin, nu_bin, t_bin, phi_bin, value))

      Q2_max = max(Q2_max, Q2_bin)
      nu_max = max(nu_max, nu_bin)
      t_max  = max(t_max, t_bin)
      phi_max = max(phi_max, phi_bin)

  acc = np.zeros((Q2_max+1, nu_max+1, t_max+1, phi_max+1))

  for Q2_bin, nu_bin, t_bin, phi_bin, value in entries:
    acc[Q2_bin, nu_bin, t_bin, phi_bin] = value

  return acc

def main():
  # Convert Johannes' txt files into arrays 
  #file_muPlus = "/afs/cern.ch/user/g/gkainth/Johannes_code/2016-dvcs-analysis/output/dvcs_acc/P04/out_acc_0.txt"
  #acc_muPlus_J = load_colleague_acc(file_muPlus)

  #file_muMinus = "/afs/cern.ch/user/g/gkainth/Johannes_code/2016-dvcs-analysis/output/dvcs_acc/P04/out_acc_1.txt"
  #acc_muMinus_J = load_colleague_acc(file_muMinus)

  # Get Johannes' arrays from the pickle
  with open("acc_J_P04.pkl", "rb") as f:
    acc_muPlus_J  = pickle.load(f)
    acc_muMinus_J = pickle.load(f)

  # Get my arrays from the pickle 
  with open("G_acc_P04.pkl", "rb") as f:
    acc_muPlus  = pickle.load(f)
    acc_muMinus = pickle.load(f)

  print("Johannes:")
  print("mu+ acceptance min:", np.min(acc_muPlus_J))
  print("mu+ acceptance max:", np.max(acc_muPlus_J))
  print("mu+ acceptance mean:", np.mean(acc_muPlus_J))
  print("mu- acceptance min:", np.min(acc_muMinus_J))
  print("mu- acceptance max:", np.max(acc_muMinus_J))
  print("mu- acceptance mean:", np.mean(acc_muMinus_J))

  print("\nGursimran:")
  print("mu+ acceptance min:", np.min(acc_muPlus))
  print("mu+ acceptance max:", np.max(acc_muPlus))
  print("mu+ acceptance mean:", np.mean(acc_muPlus))
  print("mu- acceptance min:", np.min(acc_muMinus))
  print("mu- acceptance max:", np.max(acc_muMinus))
  print("mu- acceptance mean:", np.mean(acc_muMinus))

  print("\nComparison:")
  ratio = np.divide(
    acc_muPlus,
    acc_muPlus_J,
    out=np.zeros_like(acc_muPlus),
    where=acc_muPlus_J > 0
  )

  print("ratio mean:", ratio[ratio > 0].mean())
  print("ratio min :", ratio[ratio > 0].min())
  print("ratio max :", ratio[ratio > 0].max())

  print("G count:", np.count_nonzero(acc_muPlus))
  print("J count:", np.count_nonzero(acc_muPlus_J))

  print("G count:", acc_muPlus.size)
  print("J count:", acc_muPlus_J.size)

  """ with open("acc_J_P04.pkl", "wb") as f:
    pickle.dump(acc_muPlus_J, f)
    pickle.dump(acc_muMinus_J, f)
  print("Saved to acc_J_P04.pkl.") """

if __name__ == "__main__":
  main()