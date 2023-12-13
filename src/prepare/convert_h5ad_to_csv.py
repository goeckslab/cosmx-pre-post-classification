import anndata as an
import pandas as pd
from pathlib import Path

csv_path: Path = Path("data", "csv", "0")

if not csv_path.exists():
    csv_path.mkdir(parents=True, exist_ok=True)

data: an.AnnData = an.read_h5ad("data/h5ad/BEMS265303.h5ad")
# convert h5ad to pandas dataframe
df = pd.DataFrame(data.X, index=data.obs.index, columns=data.var.index)

metadata = pd.read_excel("data/Confidential-SMMART Patient Metadata_23MAY2022.xlsx")

# iterate through data folder
for file in Path("data/h5ad").iterdir():
    if "Confidential" in file.name:
        continue

    # create new variable sampple id, extract it from the last part of the file name and drop BEMS
    sample_id = file.name.split("_")[-1].replace("BEMS", "")
    # drop the .h5ad part
    sample_id = sample_id.replace(".h5ad", "")

    data = an.read_h5ad(file)

    # convert h5ad to pandas dataframe
    df = pd.DataFrame(data.X, index=data.obs.index, columns=data.var.index)
    # find Patient ID in metadata associated with the sample id
    patient_id = metadata[metadata["Sample ID"] == int(sample_id)]["De-identified Patient ID"].values[0]

    # skip patient a and b
    if patient_id == "Patient A" or patient_id == "Patient B":
        continue

    var = df.var().mean()
    mean = df.mean().mean()


    df["Patient Id"] = patient_id
    df["Sample Id"] = sample_id

    # save dataframe to csv
    df.to_csv(Path(csv_path, f"{patient_id.replace(' ', '_').lower()}_{sample_id}.csv"), index=False)
