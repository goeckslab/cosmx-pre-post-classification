import pandas as pd
from pathlib import Path
import os
from typing import List


def load_patient(patient: str) -> {}:
    data_path = Path("data", "mapped_data")

    data_frames: dict = {}

    for root, dirs, files in os.walk(data_path):
        for file in files:

            if Path(file).suffix != ".csv":
                continue

            if patient not in file:
                continue

            print(file)

            data_frames[Path(file).stem] = pd.read_csv(Path(data_path, file))

    return data_frames


def load_files(patient_to_be_excluded: str) -> (pd.DataFrame, List):
    if not patient_to_be_excluded:
        raise ValueError("Patient to be excluded needs to be specified.")

    data_path = Path("data", "mapped_data")

    data_frames: [pd.DataFrame] = []
    loaded_files = []
    for root, dirs, files in os.walk(data_path):
        for file in files:

            if Path(file).suffix != ".csv":
                continue

            if patient_to_be_excluded in file:
                continue

            print(file)

            data_frames.append(pd.read_csv(Path(data_path, file)))
            loaded_files.append(file)

    data_frames = pd.concat(data_frames, axis=0)
    return data_frames, loaded_files


def load_all_files() -> (pd.DataFrame, List):
    data_path = Path("data", "mapped_data")

    data_frames: [pd.DataFrame] = []
    loaded_files = []
    for root, dirs, files in os.walk(data_path):
        for file in files:

            if Path(file).suffix != ".csv":
                continue

            print(file)

            data_frames.append(pd.read_csv(Path(data_path, file)))
            loaded_files.append(file)

    data_frames = pd.concat(data_frames, axis=0)
    return data_frames, loaded_files
