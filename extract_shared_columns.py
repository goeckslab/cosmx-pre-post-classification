import os
import numpy as np
import pandas as pd

ddir = 'data_v3/csv/0'

markers = []
for fh in os.listdir(ddir):
    if fh.endswith('.csv'):
        df = pd.read_csv(os.path.join(ddir, fh))
        first_protein_col_index = df.columns.get_loc('4-1BB')
        df = df.iloc[:, first_protein_col_index:]
        markers.append(set(df.columns))

marker_intersection = list(set.intersection(*markers))

print(len(marker_intersection))  # I got 107 markers
print(marker_intersection)