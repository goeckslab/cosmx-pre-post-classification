import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import numpy as np
import argparse
import scanpy as sc
from tqdm import tqdm

metadata = pd.read_excel("data/Confidential-SMMART Patient Metadata_23MAY2022.xlsx")
V2_SHARED_MARKERS = ['pRB', 'CD45', 'CK19', 'Ki67', 'aSMA', 'Ecad', 'PR', 'CK14', 'HER2', 'AR', 'CK17', 'p21',
                     'Vimentin',
                     'pERK', 'EGFR', 'ER']
V3_SHARED_MARKERS = ['VISTA', 'CD163', 'FOXP3', 'IDO1', 'STING', 'CASP9', 'ERalpha', 'pMEK1', 'panRAS', 'CD44', 'CD38',
                     'MET', 'ICAM1', 'p53', 'p44-42MAP-ERK', 'EpCAM', 'VIM', 'IgD', 'KI67', 'FABP4', 'TIM3', 'TCF-1',
                     'PDL1', 'CD4', 'E-cad', 'CD39', 'AR', 'ICOS', 'CD14', 'CD15', 'CD8', 'GATA-3', 'PDL2', 'CD3',
                     'CD66b', 'GZMA', 'GSK3b', 'FAP', 'PARP', 'AKT', 'CD19', 'PTEN', 'pTuberin', 'CD45', 'CD123',
                     'CD31', 'p38MAPK', 'INPP4B', 'CD45RA', 'CD11a', 'CCR7', 'CD138', 'CD11b', 'IL-1B', 'pJNK', 'CD56',
                     'GITR', 'CD27', 'NF-kBp65', 'LAG3', 'CD16', 'BCL2', 'pGSK3ab', 'MART-1', 'CTLA4', 'SMA', 'GZMB',
                     'BIM', 'CD68', 'p44-42MAPK-ERK', 'CD25', 'CyclinD1', 'HLA-DRA', 'CD45RO', 'OX40L', 'LAMP1',
                     'gH2AX', 'IL-18', 'BRAF', 'PLCG1', 'CD20', 'PVR', 'B7-H3', 'iNOS', '4-1BB', 'Desmin', 'CD11c',
                     'CD40', 'CD226', 'cJun', 'pPRAS40', 'PD1', 'FN1', 'MPO', 'Beta-catenin', 'CD34', 'BCLXL', 'EGFR',
                     'RSK1p90', 'Pan-AKT', 'NY-ESO-1', 'GAPDH', 'CD80', 'T-bet', 'BCL6', 'Arg', 'CD127']

if __name__ == '__main__':

    # create args parser
    parser = argparse.ArgumentParser(description='Create cell community spots.')
    parser.add_argument('--radius', "-r", type=int, help='Radius to use for spot generation', default=30,
                        choices=[15, 30, 60, 90, 120])
    parser.add_argument('--data_version', "-v", type=str, help='Version of the data to use', choices=["v1", "v2", "v3"],
                        default="v1")
    args = parser.parse_args()

    radius: int = args.radius
    data_version: str = args.data_version

    print(f"Using data version {data_version} and radius {radius}")

    save_path = Path("data_2", "csv")
    if data_version == "v1":
        save_path = Path("data", "csv")
    elif data_version == "v2":
        save_path = Path("data_2", "csv")
    elif data_version == "v3":
        save_path = Path("data_v3", "csv")

    save_path = Path(save_path, f"{radius}")

    print(f"Using save path: {save_path}")

    if not save_path.exists():
        save_path.mkdir(parents=True, exist_ok=True)

    if data_version == "v1":
        # iterate through biopsies
        biopsy_folder = Path("data", "h5ad")

        for biopsy_file in biopsy_folder.iterdir():
            if not biopsy_file.is_file():
                continue

            print(f"Processing {biopsy_file}")
            # create new variable sample id, extract it from the last part of the file name and drop BEMS
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

    elif data_version == "v2":
        biopsy_folder = Path("data_2", "csv", "0")

        for biopsy_file in biopsy_folder.iterdir():
            if not biopsy_file.is_file():
                continue

            print(f"Processing {biopsy_file}")

            df = pd.read_csv(biopsy_file, sep=",")

            proteins = df[V2_SHARED_MARKERS].copy()

            # Build the BallTree
            cell_locs = df[['X_centroid', 'Y_centroid']].values
            tree = BallTree(cell_locs)

            # Create community spot for each cell.
            communities = np.array([])
            for index, cell_loc in enumerate(cell_locs):
                # print (f"Processing cell {index} with location {cell_loc}")
                indices = tree.query_radius([cell_loc], r=radius)[0]
                # if len(indices) < 5:
                #     continue
                # print(f"Found {len(indices)} cells within distance 30 from cell {index}")
                community = proteins.iloc[indices]
                c_mean = community.mean(axis=0)

                proteins.iloc[index] = c_mean

            df = pd.DataFrame(proteins, columns=proteins.columns)
            patient_id = biopsy_file.stem
            df["Patient Id"] = patient_id

            # save biopsy
            df.to_csv(Path(save_path, f"{patient_id.replace(' ', '_').lower()}.csv"), index=False)

    elif data_version == "v3":
        biopsy_folder = Path("data_v3", "csv", "0")

        for biopsy_file in biopsy_folder.iterdir():
            if not biopsy_file.is_file():
                continue

            print(f"Processing {biopsy_file}")
            # no index column
            df = pd.read_csv(biopsy_file, sep=",", index_col=0)

            proteins = df[V3_SHARED_MARKERS].copy()

            # Build the BallTree
            cell_locs = df[['x_slide_mm', 'y_slide_mm']].values
            tree = BallTree(cell_locs)

            # Create community spot for each cell.
            communities = np.array([])
            for index, cell_loc in enumerate(cell_locs):
                # print (f"Processing cell {index} with location {cell_loc}")
                indices = tree.query_radius([cell_loc], r=radius)[0]
                # if len(indices) < 5:
                #     continue
                # print(f"Found {len(indices)} cells within distance 30 from cell {index}")
                community = proteins.iloc[indices]
                c_mean = community.mean(axis=0)

                proteins.iloc[index] = c_mean

            df = pd.DataFrame(proteins, columns=proteins.columns)
            patient_id = biopsy_file.stem
            df["Patient Id"] = patient_id

            # save biopsy
            df.to_csv(Path(save_path, f"{patient_id.replace(' ', '_').lower()}.csv"), index=False)

    else:
        raise ValueError("Unknown data version")
