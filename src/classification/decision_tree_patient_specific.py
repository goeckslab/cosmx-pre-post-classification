import pandas as pd
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
import files
import argparse
from pathlib import Path
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split


def create_results_folder(save_folder: Path) -> Path:
    experiment_id = 0
    base_path = Path(save_folder, "experiment_run")
    save_path = Path(str(base_path) + "_" + str(experiment_id))
    while Path(save_path).exists():
        save_path = Path(str(base_path) + "_" + str(experiment_id))
        experiment_id += 1

    created: bool = False
    if not save_path.exists():
        while not created:
            try:
                save_path.mkdir(parents=True)
                created = True
            except:
                experiment_id += 1
                save_path = Path(str(base_path) + "_" + str(experiment_id))

    return save_path


if __name__ == '__main__':
    save_path: Path = create_results_folder(Path("results"))
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--patient", required=True, help="Patient to be excluded", type=str)

    args = parser.parse_args()

    patient: str = args.patient
    # test_file_path: Path = args.test_file

    train_set, _ = files.load_files(patient)
    test_set = files.load_patient(patient=patient)

    train_set = train_set.drop(columns=["Patient Id", "Sample Id", "CellID", "MouseIgG1"])
    test_set = test_set.drop(columns=["Patient Id", "Sample Id", "CellID", "MouseIgG1"])
    # label encoder treatment
    le = LabelEncoder()
    train_set["Treatment"] = le.fit_transform(train_set["Treatment"])
    test_set["Treatment"] = le.fit_transform(test_set["Treatment"])

    y_train = train_set["Treatment"]
    y_test = test_set["Treatment"]

    train_set = train_set.drop(columns=["Treatment"])
    test_set = test_set.drop(columns=["Treatment"])

    min_max_scaler = MinMaxScaler(feature_range=(0, 1))
    train_set = pd.DataFrame(min_max_scaler.fit_transform(np.log10(train_set + 1)), columns=train_set.columns)
    test_set = pd.DataFrame(min_max_scaler.fit_transform(np.log10(test_set + 1)), columns=test_set.columns)

    # train decision tree classifier
    clf = DecisionTreeClassifier()
    clf.fit(train_set, y_train)
    # predict
    predictions = pd.DataFrame(data=clf.predict(test_set), columns=["Predictions"])

    # score decision tree
    accuracy = clf.score(test_set, y_test)

    data = [{"Accuracy": accuracy, "Patient": patient}]
    results = pd.DataFrame.from_records(data=data)
    results.to_csv("accuracy.csv")
    classes = pd.DataFrame()
    classes["Single Cell Id"] = y_test.index
    classes["GT"] = y_test.values
    classes["Predicted"] = predictions["Predictions"].values
    results.to_csv(Path(save_path, "accuracy.csv"), index=False)
    classes.to_csv(Path(save_path, "predictions.csv"), index=False)
