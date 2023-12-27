import os
import pandas as pd





# iterate through data_v3 folder
for root, dirs, files in os.walk("data_v3"):
    # iterate through files
    for file in files:
        # skip non csv files
        if not file.endswith(".csv"):
            continue

        file_name_parts = file.split("_")

        # rename file
        os.rename(os.path.join(root, file), os.path.join(root, f"{file_name_parts[2]}_{file_name_parts[3]}.csv").lower())

