import pandas as pd
import seaborn as sns
from pathlib import Path
import os
import matplotlib.pyplot as plt

save_path = Path("plots","clipped")
if not save_path.exists():
    save_path.mkdir(parents=True, exist_ok=True)


if __name__ == '__main__':

    data = []
    # iterate through data/outlier/2std folder
    for root, sub_directories, files in os.walk(Path("data", "outlier", "2_std")):
        for file in files:
            if 'meta' in file:
                continue
            df = pd.read_csv(Path(root, file))
            df["STD"] = 2
            # calculate percentage of clipped values

            data.append(df)

    # iterate through data/outlier/3std folder
    for root, sub_directories, files in os.walk(Path("data", "outlier", "3_std")):
        for file in files:
            if 'meta' in file:
                continue

            df = pd.read_csv(Path(root, file))
            df["STD"] = 3

            data.append(df)

    data = pd.concat(data)

    # rename Replaced to Clipped
    data = data.rename(columns={"Replaced": "Count"})

    # create bar plot for protein as x and replaced as y , STD as hue
    fig = plt.figure(figsize=(10, 5), dpi=300)
    ax = sns.boxenplot(data=data, x="Protein", y="Count", hue="STD")
    # rotate x axis laels
    plt.xticks(rotation=90)
    # strcht lgend over 2 columns
    ax.legend(ncol=2, loc="upper center")
    plt.title("Clipped Values per Protein (2 and 3 STD)")
    plt.tight_layout()
    plt.savefig(Path(save_path, "clipped_values_per_protein.png"), dpi=300)
    plt.close()


    # create bar plot for protein as x and replaced as y , STD as hue
    fig = plt.figure(figsize=(10, 5), dpi=300)
    ax = sns.boxenplot(data=data, x="Protein", y="%", hue="STD")
    # rotate x axis laels
    plt.xticks(rotation=90)
    # strcht lgend over 2 columns
    ax.legend(ncol=2, loc="upper center")
    plt.title("Clipped % per Protein (2 and 3 STD)")
    plt.tight_layout()
    plt.savefig(Path(save_path, "clipped_percentages_per_protein.png"), dpi=300)
    plt.close()
