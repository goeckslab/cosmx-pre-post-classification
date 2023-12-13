import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import numpy as np
import argparse
import scanpy as sc
from tqdm import tqdm

save_path = Path("data", "csv")
metadata = pd.read_excel("data/Confidential-SMMART Patient Metadata_23MAY2022.xlsx")

if __name__ == '__main__':

    # create args parser
    parser = argparse.ArgumentParser(description='Create cell community spots.')
    parser.add_argument('--radius', "-r", type=int, help='Radius to use for spot generation', default=30,
                        choices=[15, 30, 60, 90, 120])
    args = parser.parse_args()

    radius: int = args.radius

    save_path = Path(save_path, f"{radius}")

    if not save_path.exists():
        save_path.mkdir(parents=True, exist_ok=True)

    # iterate through biopsies
    biopsy_folder = Path("data", "h5ad")

    for biopsy_file in biopsy_folder.iterdir():
        if not biopsy_file.is_file():
            continue

        print(f"Processing {biopsy_file}")
        # create new variable sampple id, extract it from the last part of the file name and drop BEMS
        sample_id = biopsy_file.name.split("_")[-1].replace("BEMS", "")
        sample_id = sample_id.replace(".h5ad", "")
        patient_id = metadata[metadata["Sample ID"] == int(sample_id)]["De-identified Patient ID"].values[0]

        # load h5ad file
        adata = sc.read_h5ad(biopsy_file)
        proteins: pd.DataFrame = pd.DataFrame(data=adata.X, columns=adata.var_names).copy()

        spatial_information: pd.DataFrame = adata.obs[["X_centroid", "Y_centroid"]]
        spatial_information.reset_index(drop=True, inplace=True)

        # adata_copy = adata.copy()
        # Build the BallTree
        cell_locs = spatial_information[['X_centroid', 'Y_centroid']].values
        tree = BallTree(cell_locs)

        # Create community spot for each cell.
        communities = np.array([])
        for index, cell_loc in enumerate(cell_locs):
            # print (f"Processing cell {index} with location {cell_loc}")
            indices = tree.query_radius([cell_loc], r=radius)[0]
            # if len(indices) < 5:
            #     continue
            # print(f"Found {len(indices)} cells within distance 30 from cell {index}")
            community = adata.X[indices]
            c_mean = community.mean(axis=0)

            adata.X[index] = c_mean

        df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
        df["Patient Id"] = patient_id
        df["Sample Id"] = sample_id

        # save biopsy
        df.to_csv(Path(save_path, f"{patient_id.replace(' ', '_').lower()}_{sample_id}.csv"), index=False)
