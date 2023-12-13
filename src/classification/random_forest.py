import pandas as pd
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
import files
import argparse
from pathlib import Path
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

patients = ["patient_a", "patient_b", "patient_c", "patient_d"]


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
    # parser.add_argument("-p", "--patient", required=True, help="Patient to be excluded", type=str)

    args = parser.parse_args()

    # patient: str = args.patient
    # test_file_path: Path = args.test_file
    class_predictions: dict = {}
    results = []
    for patient in patients:
        train_set, _ = files.load_files(patient)
        test_sets: dict = files.load_patient(patient=patient)

        train_set = train_set.drop(columns=["Patient Id", "Sample Id", "CellID", "MouseIgG1"])

        # label encoder treatment
        le = LabelEncoder()
        train_set["Treatment"] = le.fit_transform(train_set["Treatment"])
        y_train = train_set["Treatment"]
        train_set = train_set.drop(columns=["Treatment"])
        min_max_scaler = MinMaxScaler(feature_range=(0, 1))

        train_set = pd.DataFrame(min_max_scaler.fit_transform(np.log10(train_set + 1)), columns=train_set.columns)

        # train decision tree classifier
        clf = RandomForestClassifier()
        clf.fit(train_set, y_train)

        for file_name, test_set in test_sets.items():
            test_set = test_set.drop(columns=["Patient Id", "Sample Id", "CellID", "MouseIgG1"])
            # label encoder treatment
            le = LabelEncoder()
            treatment = test_set["Treatment"].iloc[0]
            test_set["Treatment"] = le.fit_transform(test_set["Treatment"])

            y_test = test_set["Treatment"]

            test_set = test_set.drop(columns=["Treatment"])

            min_max_scaler = MinMaxScaler(feature_range=(0, 1))
            test_set = pd.DataFrame(min_max_scaler.fit_transform(np.log10(test_set + 1)), columns=test_set.columns)

            # predict
            predictions = pd.DataFrame(data=clf.predict(test_set), columns=["Predictions"])

            # score decision tree
            accuracy = clf.score(test_set, y_test)

            results.append(
                {"Accuracy": accuracy, "Patient": patient, "Model": "Random Tree", "File": file_name,
                 "Treatment": treatment})

            predicted_classes: pd.DataFrame = pd.DataFrame()
            predicted_classes["Single Cell Id"] = y_test.index
            predicted_classes["GT"] = y_test.values
            predicted_classes["Predicted"] = predictions["Predictions"].values
            class_predictions[file_name] = predicted_classes

    results = pd.DataFrame.from_records(data=results)
    results.to_csv("accuracy.csv")
    results.to_csv(Path(save_path, "accuracy.csv"), index=False)
    for file_name in class_predictions.keys():
        class_predictions[file_name].to_csv(Path(save_path, f"{file_name}_class_predictions.csv"), index=False)
