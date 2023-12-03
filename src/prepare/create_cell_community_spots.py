import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import numpy as np
import argparse
import scanpy as sc

save_path = Path("data", "csv")

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

        # load h5ad file
        adata = sc.read_h5ad(biopsy_file)
        observations = adata.obs
        # extract X_centroid and Y_centroid from observations
        spatial_information: pd.DataFrame = observations[["X_centroid", "Y_centroid"]]

        df = pd.read_csv(biopsy_file, sep="\t")
        # adata_copy = adata.copy()
        # Build the BallTree
        cell_locs = df['X_centroid', 'Y_centroid'].values
        tree = BallTree(cell_locs)

        # Create community spot for each cell.
        communities = np.array([])
        for index, cell_loc in enumerate(cell_locs):
            # print (f"Processing cell {index} with location {cell_loc}")
            indices = tree.query_radius([cell_loc], r=radius)[0]
            # if len(indices) < 5:
            #     continue
            # print(f"Found {len(indices)} cells within distance 30 from cell {index}")
            community = df[indices]
            c_mean = community.mean(axis=0)

            print(f"Community mean:")
            print(c_mean)
            df[index] = c_mean

        print(df)

        # save biopsy
        df.to_csv(Path(save_path, biopsy_file.name), sep="\t", index=False)
