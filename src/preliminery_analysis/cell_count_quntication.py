import pandas as pd
import os
from pathlib import Path

if __name__ == '__main__':
    # iterate through data/csv folder and count the number of cells per biopsy

    cell_count = {}
    for root, dirs, files in os.walk(Path("data", "csv")):
        for file in files:
            if not file.endswith(".csv"):
                continue
            if Path(root).name != "csv":
                continue

            if "csv_data" in file:
                continue

            df = pd.read_csv(Path(root, file))
            cell_count[file] = len(df)
            print(f"Found {len(df)} cells in {file}")

    cell_count = pd.DataFrame.from_dict(cell_count, orient="index", columns=["Cell Count"])
    cell_count.reset_index(inplace=True)
    cell_count.rename(columns={"index": "Biopsy"}, inplace=True)
    # print total cell count
    print(f"Total cell count: {cell_count['Cell Count'].sum()}")