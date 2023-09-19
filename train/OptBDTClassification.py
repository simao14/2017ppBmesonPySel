import numpy as np
import optuna
import pandas as pd
import uproot
from sklearn.model_selection import train_test_split
import xgboost as xgb
import os
from sklearn import metrics
import sys
sys.path.insert(0, '../')
import utils

ptmin = 3
ptmax = 60
iter = 20
varstage = [0,2,4,7,8,9,11,15] #0,2,4,7,8,11

def objective(trial):

    fileS = uproot.open("/lstore/cms/simao/sample/BPMC_3_60_small2.root")
    fileB = uproot.open("/lstore/cms/simao/sample/BPData_3_60_small2.root")

    signal,background = utils.prepdata(fileS,fileB,ptmin,ptmax)
    stage=utils.varset(varstage)

    signal_var=signal[stage]
    background_var=background[stage]

    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    
    train_x, valid_x, train_y, valid_y = train_test_split(x.to_numpy(), y.to_numpy(), test_size=0.2)
    dtrain = xgb.DMatrix(train_x, label=train_y)
    dvalid = xgb.DMatrix(valid_x, label=valid_y)
    dall = xgb.DMatrix(x.to_numpy(), label=y.to_numpy())

    param = {
        "verbosity": 0,
        "objective": "binary:logistic",
        "eval_metric": "auc",
        "booster": trial.suggest_categorical("booster", ["gbtree", "dart"]),
        "lambda": trial.suggest_float("lambda", 1e-8, 1.0, log=True),
        "alpha": trial.suggest_float("alpha", 1e-8, 1.0, log=True),
    }

    if param["booster"] == "gbtree" or param["booster"] == "dart":
        param["max_depth"] = trial.suggest_int("max_depth", 1, 20)
        param["eta"] = trial.suggest_float("eta", 1e-8, 1.0, log=True)
        param["gamma"] = trial.suggest_float("gamma", 1e-8, 1.0, log=True)
        param["grow_policy"] = trial.suggest_categorical("grow_policy", ["depthwise", "lossguide"])
    if param["booster"] == "dart":
        param["sample_type"] = trial.suggest_categorical("sample_type", ["uniform", "weighted"])
        param["normalize_type"] = trial.suggest_categorical("normalize_type", ["tree", "forest"])
        param["rate_drop"] = trial.suggest_float("rate_drop", 1e-8, 1.0, log=True)
        param["skip_drop"] = trial.suggest_float("skip_drop", 1e-8, 1.0, log=True)

    # Add a callback for pruning.
    #pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "validation-auc")
    #bst = xgb.train(param, dtrain, iter, evals=[(dvalid, "validation")], callbacks=[pruning_callback])
    #yscore = bst.predict(dvalid)
    #nn_fpr, nn_tpr, _ = metrics.roc_curve(valid_y, yscore)
    #auc=metrics.auc(nn_tpr,1-nn_fpr)

    pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "test-auc")
    history = xgb.cv(param, dall, iter, callbacks=[pruning_callback])

    auc = history["test-auc-mean"].values[-1]
    return auc


if __name__ == "__main__":
    if not os.path.exists("../results/models"):
        os.makedirs("../results/models")
    study_name="BDToptim.study"
    storage_name = "sqlite:///../results/models/{}_{}_{}.db".format(study_name,ptmin,ptmax)
    study = optuna.create_study(study_name=study_name, storage=storage_name,
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=5), direction="maximize"
    )
    study.optimize(objective, n_trials=100)
   
    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))

    print("Best trial:")
    trial = study.best_trial

    print("  Value: ", trial.value)

    print("  Params: ")
    for key, value in trial.params.items():
        print("    {}: {}".format(key, value))