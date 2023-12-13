import pandas as pd
import os, argparse
from pathlib import Path
from mappings import BIOPSY_MAPPINGS

save_path = Path("data", "mapped_data")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create cell community spots.')
    parser.add_argument('--radius', "-r", type=int, help='Radius to use for spot generation', default=30,
                        choices=[0, 15, 30, 60, 90, 120])
    args = parser.parse_args()

    radius: int = args.radius

    save_path = Path(save_path, f"{radius}")

    if not save_path.exists():
        save_path.mkdir(parents=True)

    data_path = Path("data", "csv", str(radius))
    for root, dirs, files in os.walk(data_path):
        for file in files:
            print(f"Preparing file {file}")
            patient_id = Path(file).stem.split('_')[-1]

            if patient_id in BIOPSY_MAPPINGS:
                # load dataframe
                df = pd.read_csv(Path(root, file))

                df["Treatment"] = BIOPSY_MAPPINGS[patient_id]

                assert "Treatment" in df.columns, f"Treatment column is missing for dataframe of patient {patient_id}"
                df.to_csv(Path(save_path, f"{file}"), index=False)
