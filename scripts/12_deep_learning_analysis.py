import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from keras.layers import Dense
from keras.models import Sequential
from keras.callbacks import EarlyStopping
# Import test split
from sklearn.model_selection import train_test_split
from keras import backend as K
import tensorflow as tf

pileup = bool(sys.argv[1])

if pileup == False:
    col_names = ["read_ID"] + [f"pred{i}" for i in range(1,22)] + ["type"]
elif pileup == True:
    col_names = [f"pred{i}" for i in range(1, 101)]
    col_names.append("mapq")
    col_names.append("type")


# Read datasets and skip first 5 rows that correspond to GEM output.
dataset_0 = pd.read_csv(f"""./data/machine_learning/predictors/common.{"nt_" if pileup == False else ""}predictors""", sep = "\t")
dataset_1 = pd.read_csv(f"""./data/machine_learning/predictors/novel.{"nt_" if pileup == False else ""}predictors""", sep = "\t")
# Add a column for the type of read: 0 or 1.
dataset_0["type"] = 0
dataset_1["type"] = 1
# Change column names for both datasets.
dataset_0.columns = col_names
dataset_1.columns = col_names
# Merge datasets.
master_dataset = pd.concat([dataset_0, dataset_1])
# Separate dataset in features and target variables.
feature_cols = col_names[1:-1]
target_col = "type"
x = master_dataset[feature_cols]
y = master_dataset[target_col]
# Split data: 0.75 training, 0.25 test.
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.25, random_state = 0)

model = Sequential()
model.add(Dense(101, activation = "relu", input_shape = (x_test.shape[1],)))
for i in range(50):
    model.add(Dense(500, activation = "relu"))
model.add(Dense(1, activation = "softmax")) # Best for categorical outputs.

model.compile(loss = "binary_crossentropy", optimizer = "adam", metrics = ["accuracy"])

model.fit(x_train, y_train, epochs = 20)