import os
import uproot
import optuna
from   optuna.trial import TrialState
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import metrics

import sys
sys.path.insert(0, '../')
import utils

DEVICE = torch.device("cpu")
BATCHSIZE = 256
CLASSES = 2
EPOCHS = 20

ptmin=3
ptmax=60

varstage=[0,2,4,7,8,9,11,15]

def define_model(trial):
    # We optimize the number of layers, hidden units and dropout ratio in each layer.
    n_layers = trial.suggest_int("n_layers", 1, 5)
    layers = []

    in_features = len(varstage)
    for i in range(n_layers):
        out_features = trial.suggest_int("n_units_l{}".format(i), 128, 256)
        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = trial.suggest_float("dropout_l{}".format(i), 0.0, 0.5)
        layers.append(nn.Dropout(p))

        in_features = out_features
    layers.append(nn.Linear(in_features, CLASSES))

    return nn.Sequential(*layers)


def objective(trial):
    # Generate the model.
    model = define_model(trial).to(DEVICE)
    
    # Generate the optimizers.
    optimizer_name = trial.suggest_categorical("optimizer", ["Adam", "RMSprop", "SGD"])
    lr = trial.suggest_float("lr", 1e-6, 1e-2, log=True)
    wd = trial.suggest_float("weight_decay", 1e-4, 0.05)
    optimizer = getattr(optim, optimizer_name)(model.parameters(), lr=lr, weight_decay=wd)

   
    fileS = uproot.open("/lstore/cms/simao/sample/BPMC_3_60_small2.root")
    fileB = uproot.open("/lstore/cms/simao/sample/BPData_3_60_small2.root")

    signal,background = utils.prepdata(fileS,fileB,ptmin,ptmax)
    stage=utils.varset(varstage)

    signal_var=signal[stage]
    background_var=background[stage]

    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    dataset = utils.ClassificationDataset(x.values,y.values)
    train_set,val_set =torch.utils.data.random_split(dataset, [0.8,0.2])
    train_loader = torch.utils.data.DataLoader(train_set, batch_size=BATCHSIZE, shuffle=True)
    valid_loader = torch.utils.data.DataLoader(val_set, batch_size=BATCHSIZE, shuffle=True)
    dev_X, dev_y = val_set[:][0], val_set[:][1]
    # Training of the model.
    for epoch in range(EPOCHS):
        model.train()
        for data, target in train_loader:
            
            data, target = data.view(data.size(0), -1).to(DEVICE), target.to(DEVICE)

            optimizer.zero_grad()
            output = model(data)
            loss = F.cross_entropy(output, target)
            loss.backward()
            optimizer.step()

        # Validation of the model.
        model.eval()
        correct = 0
        with torch.no_grad():
            for data, target in valid_loader:
                
                data, target = data.view(data.size(0), -1).to(DEVICE), target.to(DEVICE)
                output = model(data)
                # Get the index of the max log-probability.
                pred = output.argmax(dim=1, keepdim=True)
                correct += pred.eq(target.view_as(pred)).sum().item()

            #data_x=pd.DataFrame(dev_X.numpy())
            #data_x.columns=stage
            #data_x["label"]=dev_y
            #data_x["DNN_output"]=torch.nn.functional.softmax(model(dev_X),dim=1)[:,1].detach().numpy()
            #data_x["scale_weight"]=np.ones(len(dev_y))
            #sig,_,_,_=utils.significancecalc(data_x,100,"DNN_output")

            yscore = torch.nn.functional.softmax(model(dev_X),dim=1)[:,1]
            nn_fpr, nn_tpr, nn_thresholds = metrics.roc_curve(dev_y.detach().numpy(), yscore.detach().numpy())
            auc=metrics.auc(nn_tpr,1-nn_fpr)

        #accuracy = correct / len(valid_loader.dataset)

        trial.report(auc, epoch)

        # Handle pruning based on the intermediate value.
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return auc


if __name__ == "__main__":
    if not os.path.exists("../results/models"):
        os.makedirs("../results/models")
    study_name="NNoptim.study"
    storage_name = "sqlite:///../results/models/{}_{}_{}.db".format(study_name,ptmin,ptmax)
    study = optuna.create_study(study_name=study_name, storage=storage_name,direction="maximize",load_if_exists=True)
    study.optimize(objective, n_trials=100)
    pruned_trials = study.get_trials(deepcopy=False, states=[TrialState.PRUNED])
    complete_trials = study.get_trials(deepcopy=False, states=[TrialState.COMPLETE])

    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    print("Best trial:")
    trial = study.best_trial

    print("  Value: ", trial.value)

    print("  Params: ")
    for key, value in trial.params.items():
        print("    {}: {}".format(key, value))