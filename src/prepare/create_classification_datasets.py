import pandas as pd
import os
from pathlib import Path
from mappings import BIOPSY_PRE_OR_POST

save_path = Path("data", "mapped_data")

if __name__ == '__main__':
    if not save_path.exists():
        save_path.mkdir(parents=True)

    data_path = Path("data", "csv")
    for root, dirs, files in os.walk(data_path):
        for file in files:
            patient_id = Path(file).stem.split('_')[-1]

            if patient_id in BIOPSY_PRE_OR_POST:
                print(f"Loading file {file}...")
                # load dataframe
                df = pd.read_csv(Path(root, file))

                df["Treatment"] = BIOPSY_PRE_OR_POST[patient_id]

                assert "Treatment" in df.columns, f"Treatment column is missing for dataframe of patient {patient_id}"
                print(f"Saving file {file}...")
                df.to_csv(Path(save_path, f"{file}"), index=False)
