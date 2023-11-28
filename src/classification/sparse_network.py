from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.models import Sequential
from pathlib import Path
import argparse
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
import files
import numpy as np

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
    class_predictions: dict = {}
    results = []
    for patient in patients:
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

        # Create a Sequential model
        model = Sequential([
            Dense(128, activation='relu', input_shape=(len(train_set.columns),)),
            Dropout(0.3),
            Dense(64, activation='relu'),
            Dropout(0.3),  # Dropout layer after the first Dense layer
            Dense(2, activation='softmax')  # Output layer with softmax activation for multi-class classification
        ])

        # Compile the model
        model.compile(optimizer='adam',
                      loss='sparse_categorical_crossentropy',  # Using sparse categorical crossentropy
                      metrics=['accuracy'])

        # Summary of the model
        model.summary()

        # Assuming X_train and y_train are your training data and labels
        history = model.fit(train_set, y_train, epochs=10, validation_split=0.2)

        evaluation = model.evaluate(test_set, y_test)

        results.append({"Accuracy": evaluation[1], "Patient": patient, "Model": "Sparse Neural Network"})

        # predict
        predictions = model.predict(test_set)
        predictions = pd.DataFrame(data=np.argmax(predictions, axis=1), columns=["Predictions"])

        predicted_classes: pd.DataFrame = pd.DataFrame()
        predicted_classes["Single Cell Id"] = y_test.index
        predicted_classes["GT"] = y_test.values
        predicted_classes["Predicted"] = predictions["Predictions"].values
        class_predictions[patient] = predicted_classes

    results = pd.DataFrame.from_records(data=results)
    results.to_csv("accuracy.csv")
    results.to_csv(Path(save_path, "accuracy.csv"), index=False)
    for patient in class_predictions.keys():
        class_predictions[patient].to_csv(Path(save_path, f"{patient}_class_predictions.csv"), index=False)
