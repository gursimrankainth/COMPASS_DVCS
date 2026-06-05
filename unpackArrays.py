import pickle
import numpy as np

# ***************************
def print_data(period, array_name, data_type):
  """
  data_type: "arrays", "acc", "vars", or "accVar"
  """

  filename = f"dvcs_{data_type}_{period}.pkl"

  with open(filename, "rb") as f:
    data = pickle.load(f)

  if array_name not in data:
    print("Available keys:", data.keys())
    return

  arr = data[array_name]

  for i in range(arr.shape[0]):  # t
    for j in range(arr.shape[1]):  # Q2
      for k in range(arr.shape[2]):  # nu
        for l in range(arr.shape[3]):  # phi

          # reverse t-bin ordering
          rev_i = (arr.shape[0] - 1) - i
          value = arr[i, j, k, l]

          print(value)

# ***************************
#"real_muPlus", "real_muMinus"
#"BH_muPlus", "BH_muMinus"
#"lepPi0_muPlus", "lepPi0_muMinus"
#"hepPi0_muPlus", "hepPi0_muMinus"

#print_data("P04", "BH_muPlus", "arrays")
#print_data("P04", "BH_muPlus", "vars")
#print_data("P04", "gen_muPlus", "acc")
print_data("P04", "rec_muPlus", "accVar")

