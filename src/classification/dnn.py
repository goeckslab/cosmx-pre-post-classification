import pandas as pd
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
import files
import argparse
from pathlib import Path
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.models import Model

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("-p", "--patient", required=True, help="Patient to be excluded", type=str)
    # parser.add_argument('-f', "--test_file", required=True, type=Path)

    args = parser.parse_args()

    # patient_to_be_excluded: str = args.patient
    # test_file_path: Path = args.test_file

    df, _ = files.load_all_files()
    df = df.drop(columns=["Patient Id", "Sample Id", "CellID", "MouseIgG1"])
    # label encoder treatment
    le = LabelEncoder()
    df["Treatment"] = le.fit_transform(df["Treatment"])

    y = df["Treatment"]
    X = df.drop(columns=["Treatment"])

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

    min_max_scaler = MinMaxScaler(feature_range=(0, 1))
    X_train = pd.DataFrame(min_max_scaler.fit_transform(np.log10(X_train + 1)), columns=X_train.columns)
    X_test = pd.DataFrame(min_max_scaler.fit_transform(np.log10(X_test + 1)), columns=X_test.columns)

    # build model

    model: Model = Model()

