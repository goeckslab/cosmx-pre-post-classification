import sys

import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

plot_save_path = Path("plots", "variance", "clipped", "3_std")
data_save_path = Path("data", "outlier",  "3_std")
if not plot_save_path.exists():
    plot_save_path.mkdir(parents=True, exist_ok=True)

if not data_save_path.exists():
    data_save_path.mkdir(parents=True, exist_ok=True)

columns = ['41BB', 'B7H3', 'Bcl2', 'BetaCatenin', 'CCR7', 'CD11b', 'CD11c',
           'CD127', 'CD14', 'CD19', 'CD20', 'CD27', 'CD31', 'CD34', 'CD39', 'CD3',
           'CD40', 'CD45', 'CD45RA', 'CD4', 'CD56', 'CD68', 'CD8', 'CTLA4', 'EGFR',
           'EpCAM', 'FABP4', 'Fibronectin', 'GITR', 'GZMA', 'GZMB', 'HER2',
           'Histone', 'HLADR', 'ICOS', 'IDO1', 'IgD', 'IL18', 'IL1b', 'iNOS',
           'Ki67', 'LAG3', 'MouseIgG1', 'NFkBp65', 'panRAS', 'PD1', 'PDL1', 'PDL2',
           'SMA', 'STING', 'Syndecan1', 'Tbet', 'TCF1TCF7', 'TIM3', 'Vimentin', "Patient Id", "Sample Id"]

patients = {}

meta_data = []

# iterate over csv files in data/csv
for file in Path("data/csv").iterdir():
    try:
        # load file
        df = pd.read_csv(file)[columns]

        patient_id = df["Patient Id"].unique()[0]
        sample_id = df["Sample Id"].unique()[0]
        df = df.drop(columns=["Patient Id", "Sample Id"])

        # convert all columns to int
        df = df.astype(float)

        # calculate how many values are greater than 3 std
        clipped = np.sum(np.abs(df - df.mean()) > (3 * df.std()))
        # print first columns of clipped

        clipped_per_columns = pd.DataFrame(clipped)
        # rename 0 to Replaced
        clipped_per_columns = clipped_per_columns.rename(columns={0: "Replaced"})
        # reset index
        clipped_per_columns = clipped_per_columns.reset_index()
        # rename index to Protein
        clipped_per_columns = clipped_per_columns.rename(columns={"index": "Protein"})
        clipped_per_columns["%"] = round((clipped_per_columns["Replaced"] / df.shape[0]) * 100, 2)
        # save clipped per columns
        clipped_per_columns.to_csv(Path(data_save_path, f"{sample_id}_clipped.csv"), index=False)

        try:
            # clip outliers that are greater than 3 std
            df = df.clip(lower=df.mean() - (3 * df.std()), upper=df.mean() + (3 * df.std()), axis=1)
        except BaseException as ex:
            print(ex)
            input()

        meta_data.append({
            "Sample Id": sample_id,
            "Clipped Values": clipped.sum(),
            "Total Values": df.shape[0] * df.shape[1],
            "% Clipped": round((clipped.sum() / (df.shape[0] * df.shape[1])) * 100, 2)
        })

        # scale data to 0-1
        df = (df - df.min()) / (df.max() - df.min())

        if patient_id in patients:
            patients[patient_id].append(df)
        else:
            patients[patient_id] = [df]

    except KeyboardInterrupt:
        sys.exit()

    except BaseException as ex:
        print(ex)
        print(file)

meta_data = pd.DataFrame(meta_data)
meta_data.to_csv(Path(data_save_path, "meta_data.csv"), index=False)

for patient in patients:
    # merge dfs together
    df = pd.concat(patients[patient])
    # create violin plot
    fig = plt.figure(figsize=(10, 5), dpi=300)
    sns.violinplot(data=df)
    # rotate x axis labels by 90 degree
    plt.xticks(rotation=90)
    plt.title(f"Variance of {patient}")
    plt.tight_layout()
    # save plot
    plt.savefig(Path(plot_save_path, f"{patient}.png"), dpi=300)
    plt.close()
