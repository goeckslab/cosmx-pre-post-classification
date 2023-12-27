import os
import pandas as pd



dfs = []


for root, dirs, files in os.walk("data_v3/csv/0"):
    # iterate through files
    for file in files:
        # skip non csv files
        if not file.endswith(".csv"):
            continue

        # load all files into a list of dataframes
        df = pd.read_csv(os.path.join(root, file))
        dfs.append(df)


# find the shared columsn between all dataframes
shared_columns = list(set.intersection(*map(set, dfs)))
print(shared_columns)
print(len(shared_columns))
