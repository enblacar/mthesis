import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics

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
# Initialize the model.
logreg = LogisticRegression()
# Fit the model with the data.
logreg.fit(x_train, y_train)
# Make a prediction on the test set.
y_pred = logreg.predict(x_test)
# Make a confussion matrix to evaluate the performance.
"""A confussion matrix is a table that is used to evaluate the performance of a classification model.
The fundamental of a confussion matrix is that the number of correct and incorrect predictions are summed up class-wise."""

cnf_matrix = metrics.confusion_matrix(y_test, y_pred)

## Visualize the confussion matrix in a heatmap.
class_names = [0, 1] # Name of the classes.
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
sns.heatmap(pd.DataFrame(cnf_matrix), annot = True, cmap = "YlGnBu", fmt = "g")
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title("Confussion matrix", y = 1.1)
plt.ylabel("Actual label")
plt.xlabel("Predicted label")
plt.savefig(f"""./results/{"nt_" if pileup == False else "pileup_"}cnf_matrix.png""")
plt.show()

with open(f"""./results/logreg_{"nt_" if pileup == False else "pileup_"}metrics.txt""", "wt") as out:
    out.write(f"""Accuracy: {metrics.accuracy_score(y_test, y_pred)}\n""")
    out.write(f"""Precision: {metrics.precision_score(y_test, y_pred)}\n""")
    out.write(f"""Recall: {metrics.recall_score(y_test, y_pred)}\n""")

# ROC curve
"""Receiver Operating Characteristic (ROC) curve is a plot of the true positive rate against
the false positive rate. It shows the tradeoff between sensitivity and specificity."""
y_pred_proba = logreg.predict_proba(x_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)
plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
plt.legend(loc=4)
plt.savefig(f"./results/{"nt_" if pileup == False else "pileup_"}roc_curve.png")
plt.show() # AUC score of 1 is a perfect classifier. 0.5 is a worthless classifier.