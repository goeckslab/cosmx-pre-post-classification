import pandas as pd
import os, argparse
from pathlib import Path
from mappings import BIOPSY_MAPPINGS

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create cell community spots.')
    parser.add_argument("--data_version", "-v", type=str, help="Version of the data to use", choices=["v1", "v2"])
    parser.add_argument('--radius', "-r", type=int, help='Radius to use for spot generation', default=30,
                        choices=[0, 15, 30, 60, 90, 120])
    args = parser.parse_args()

    radius: int = args.radius
    data_version: str = args.data_version

    if data_version == "v1":
        save_path = Path("data", "mapped_data")
    else:
        save_path = Path("data_2", "mapped_data")

    save_path = Path(save_path, f"{radius}")

    if not save_path.exists():
        save_path.mkdir(parents=True)

    if data_version == "v1":
        data_path = Path("data", "csv", str(radius))
    else:
        data_path = Path("data_2", "csv", str(radius))

    print(f"Using path: {data_path}")

    for root, dirs, files in os.walk(data_path):
        for file in files:
            print(f"Preparing file {file}")
            if data_version == "v1":
                patient_id = Path(file).stem.split('_')[-1]
            else:
                patient_id = Path(file).stem

            if patient_id in BIOPSY_MAPPINGS:
                # load dataframe
                df = pd.read_csv(Path(root, file))

                df["Treatment"] = BIOPSY_MAPPINGS[patient_id]

                assert "Treatment" in df.columns, f"Treatment column is missing for dataframe of patient {patient_id}"
                df.to_csv(Path(save_path, f"{file}"), index=False)
