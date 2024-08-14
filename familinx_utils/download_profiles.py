import csv
import numpy as np
import pandas as pd
from tqdm import tqdm

id_var_ids = [0,1]
birth_var_ids = [14, 19, 20, 21, 23, 24]
death_var_ids = [29, 34, 35, 36, 38, 39, 74]

with open('../familinx/profiles-anon.txt') as csvfile:
    data = []
    complete_profiles = 0
    for i, row in enumerate(pbar := tqdm(csv.reader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE))):
        data_row = []
        for j, x in enumerate(row):
            # Read header.
            if i == 0:
                data_row.append(x)
            # Skip if missing basic info.
            elif (j in id_var_ids) and (x == '*'):
                data_row = []
                break
            # Skip if missing birth data.
            elif (j in birth_var_ids) and (x == '*'):
                data_row = []
                break
            # Skip if missing death data and unlikely to be alive.
            elif (j in death_var_ids) and (int(row[2]) < 2015 - 150) and (x == '*'):
                data_row = []
                break
            # Include if passes all checks.
            else:
                data_row.append(np.nan if x == '*' else x)
        if len(data_row) != 0:
            data.append(data_row)
            complete_profiles += 1
        pbar.set_description(f"Processing {i+1} profiles ({np.round(min(i+1,68588512)/68588512 * 100, 2)}%) | Complete profiles so far: {complete_profiles} ({np.round(complete_profiles/(i+1),2 * 100)}%)")
print(len(data))
df = pd.DataFrame(data=data[1:], columns=data[0])
df.to_csv('../familinx/complete-v2-profiles-anon.txt', sep='\t')