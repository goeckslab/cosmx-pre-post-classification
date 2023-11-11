import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

save_path = Path("plots", "variance", "3_std")
if not save_path.exists():
    save_path.mkdir(parents=True, exist_ok=True)

columns = ['41BB', 'B7H3', 'Bcl2', 'BetaCatenin', 'CCR7', 'CD11b', 'CD11c',
           'CD127', 'CD14', 'CD19', 'CD20', 'CD27', 'CD31', 'CD34', 'CD39', 'CD3',
           'CD40', 'CD45', 'CD45RA', 'CD4', 'CD56', 'CD68', 'CD8', 'CTLA4', 'EGFR',
           'EpCAM', 'FABP4', 'Fibronectin', 'GITR', 'GZMA', 'GZMB', 'HER2',
           'Histone', 'HLADR', 'ICOS', 'IDO1', 'IgD', 'IL18', 'IL1b', 'iNOS',
           'Ki67', 'LAG3', 'MouseIgG1', 'NFkBp65', 'panRAS', 'PD1', 'PDL1', 'PDL2',
           'SMA', 'STING', 'Syndecan1', 'Tbet', 'TCF1TCF7', 'TIM3', 'Vimentin', "Patient Id"]

# columns = ['Ki67', 'LAG3', 'MouseIgG1', 'NFkBp65', 'panRAS', 'PD1', 'PDL1', 'PDL2',
#           'SMA', 'STING', 'Syndecan1', 'Tbet', 'TCF1TCF7', 'TIM3', 'Vimentin', "Patient Id"]

patients = {}

# iterate over csv files in data/csv
for file in Path("data/csv").iterdir():
    try:
        # load file
        df = pd.read_csv(file)[columns]

        patient_id = df["Patient Id"].unique()[0]
        df = df.drop(columns=["Patient Id"])

        # convert all columns to int
        df = df.astype(float)

        # remove outliers which are greather than 3 std
        df = df[np.abs(df - df.mean()) <= (3 * df.std())]
        # scale data to 0-1
        df = (df - df.min()) / (df.max() - df.min())

        if patient_id in patients:
            patients[patient_id].append(df)
        else:
            patients[patient_id] = [df]
    except BaseException as ex:
        print(ex)
        print(file)

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
    plt.savefig(Path(save_path, f"{patient}.png"), dpi=300)
    plt.close()
