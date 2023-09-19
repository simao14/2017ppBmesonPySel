
import argparse
import torch
import random
import uproot
import torch.nn as nn
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os.path
from sklearn import metrics

import sys
sys.path.insert(0, '../')
import utils
import utils_factors

class FeedforwardNetwork(nn.Module):
    def __init__(
            self, n_classes, n_features, hidden_size, layers,
            activation_type, dropout, **kwargs):
        
        super().__init__()

        self.first_layer = nn.Linear(n_features, hidden_size[0])
        self.hidden_layers = nn.ModuleList([nn.Linear(hidden_size[i], hidden_size[i+1]) for i in range(layers-1)]) #creates a list which holds modules, letting us automatically create different nn.Linears for the hidden layers.
        self.output_layer = nn.Linear(hidden_size[-1], n_classes)
        if activation_type == "relu": self.activation = nn.ReLU()
        else: self.activation = nn.Tanh()
        self.order = nn.Sequential(self.first_layer,self.activation,nn.Dropout(p=dropout[0]))   
        for k in range(layers-1):                                           
            self.order.append(self.hidden_layers[k])
            self.order.append(self.activation)
            self.order.append(nn.Dropout(p=dropout[k+1]))
        self.order.append(self.output_layer)                                

        
    def forward(self, x, **kwargs):
        
        x = self.order(x)                             
        return x


def train_batch(X, y, model, optimizer, criterion, **kwargs):

    model.train()
    optimizer.zero_grad()    # reset the grad value
    y2 = model(X)            # predicted scores
    loss = criterion(y2,y)   # needs to be (pred_scores,gold labels), y2 has an extra dimension due to tracking of grad
    loss.backward()          # calculates the grads
    optimizer.step()         # updates weight 
    return loss.item()       # only want the value, giving everything occupies too much memory
    


def predict(model, X):

    model.eval()
    scores = model(X)  # (n_examples x n_classes)
    predicted_labels = scores.argmax(dim=-1)  # (n_examples)
    return predicted_labels


def evaluate(model, X, y):
   
    
    y_hat = predict(model, X)
    n_correct = (y == y_hat).sum().item()
    n_possible = float(y.shape[0])
    
    return n_correct / n_possible


def plot(epochs, plottable, xlabel="Epochs", ylabel='', name='',label=None):
    plt.clf()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(epochs, plottable, label=label)
    plt.legend()
    plt.grid(visible=True)
    plt.savefig('../results/train/%s.pdf' % (name), bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ptmin', type=float)
    parser.add_argument('ptmax', type=float)
    parser.add_argument('-epochs', default=20, type=int)
    parser.add_argument('-batch_size', default=256, type=int)
    parser.add_argument('-hidden_size', type=int,nargs='+', default=[159,164,182,205])
    parser.add_argument('-layers', type=int, default=4)
    parser.add_argument('-l2_decay', type=float, default=0.0001239992504598163)
    parser.add_argument('-learning_rate', type=float, default= 0.001267306794677091)
    parser.add_argument('-dropout', type=float, nargs='+', default=[0.1380354038039047, 0.10421511990259287, 0.3093548916680354, 0.11584480748125763])
    parser.add_argument('-stages', type=int, nargs='+', default=[0,2,4,7,8,9,11,15])
    parser.add_argument('-activation',
                        choices=['tanh', 'relu'], default='relu')
    parser.add_argument('-optimizer',
                        choices=['sgd', 'adam',"rms"], default='adam')
    opt = parser.parse_args()
    
    print("files")
    fileS = uproot.open("/lstore/cms/simao/sample/BPMC_3_60_small2.root")
    fileB = uproot.open("/lstore/cms/simao/sample/BPData_3_60_small2.root")
    #fileS=uproot.open("~/Desktop/UNI/LIP/mnt/data/BPMC_3_60.root")
    #fileB=uproot.open("~/Desktop/UNI/LIP/mnt/data/BPData_3_60.root")
    
    signal,background = utils.prepdata(fileS,fileB,opt.ptmin,opt.ptmax)
    
    stage=utils.varset(opt.stages)
    signal_var=signal[stage]
    background_var=background[stage]
    
    fullstages=[0,2,4,7,8,9,11,12,14,15]
    fullstagename=utils.varset(fullstages)
    
    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    
    

    dataset = utils.ClassificationDataset(x.values,y.values)
    train_set,val_set,test_set =torch.utils.data.random_split(dataset, [0.8,0.1,0.1])
    train_dataloader = torch.utils.data.DataLoader(train_set, batch_size=opt.batch_size, shuffle=True)
    train_X, train_y = train_set[:][0], train_set[:][1]
    dev_X, dev_y = val_set[:][0], val_set[:][1]
    test_X, test_y = test_set[:][0], test_set[:][1]

    n_classes = torch.unique(dataset.y).shape[0]  #2
    n_feats = dataset.X.shape[1]  # len(stages)

    # initialize the model

    model = FeedforwardNetwork(
        n_classes,
        n_feats,
        opt.hidden_size,
        opt.layers,
        opt.activation,
        opt.dropout
    )

    # get an optimizer
    optims = {"adam": torch.optim.Adam, "sgd": torch.optim.SGD, "rms": torch.optim.RMSprop}

    optim_cls = optims[opt.optimizer]
    optimizer = optim_cls(
        model.parameters(),
        lr=opt.learning_rate,
        weight_decay=opt.l2_decay)

    # get a loss criterion
    criterion = nn.CrossEntropyLoss()

    # training loop
    epochs = torch.arange(1, opt.epochs + 1)
    train_mean_losses = []
    valid_accs = []
    train_losses = []
    train_accs = []
    valid_losses = []
    for ii in epochs:

        print('Training epoch {}'.format(ii))
        for X_batch, y_batch in train_dataloader:
            loss = train_batch(
                X_batch, y_batch, model, optimizer, criterion)
            train_losses.append(loss)

        mean_loss = torch.tensor(train_losses).mean().item()
        print('Training loss: %.4f' % (mean_loss))

        train_mean_losses.append(mean_loss)
        train_accs.append(evaluate(model, train_X, train_y))

        valid_losses.append(criterion(model(dev_X),dev_y).item())
        valid_accs.append(evaluate(model, dev_X, dev_y))
        
        print('Valid acc: %.4f' % (valid_accs[-1]))
    # plot

    if not os.path.exists("../results/models"):
        os.makedirs("../results/models")
    
    save_model = torch.jit.script(model)
    save_model.save(f"../results/models/Py_{opt.ptmin}_{opt.ptmax}.pt")

    stagelist=utils.replacespecial(opt.stages)
    config = "NN-{}-{}-{}".format(opt.ptmin, opt.ptmax,stagelist)

    if not os.path.exists("../results/train/%s" % config):
        os.makedirs("../results/train/%s" % config)
   
    plt.clf()
    plt.xlabel("Epochs")
    plt.ylabel('Loss')
    plt.plot(epochs, train_mean_losses, label="training loss")
    plt.plot(epochs, valid_losses, label="validation loss")
    plt.legend()
    plt.savefig('../results/train/%s/NN-loss.pdf' % (config), bbox_inches='tight')

    plt.clf()
    plt.xlabel("Epochs")
    plt.ylabel('Accuracy')
    plt.plot(epochs, train_accs, label="training accuracy")
    plt.plot(epochs, valid_accs, label="validation accuracy")
    plt.legend()
    plt.savefig('../results/train/%s/NN-accuracy.pdf' % (config), bbox_inches='tight')

    utils.varplot(stage,signal_var,background_var,config,"train")

    yscore = torch.nn.functional.softmax(model(test_X),dim=1)[:,1]
    nn_fpr, nn_tpr, nn_thresholds = metrics.roc_curve(test_y.detach().numpy(), yscore.detach().numpy())
    auc=metrics.auc(nn_tpr,1-nn_fpr)
    plot(nn_tpr,1-nn_fpr,ylabel="background efficiency", xlabel="signal efficiency",name='{}/ROC-curve'.format(config),label=f'auc:{auc:.3f};')

    data_x=pd.DataFrame(test_X.numpy())
    data_x.columns=stage
    data_x["label"]=test_y
    data_x["DNN_output"]=yscore.detach().numpy()
    
    data_x["scale_weight"]=np.ones(len(test_y))

    variable="DNN_output"
    fs,fb = utils_factors.get_factors(fileS,uproot.open("/user/s/smcosta/data/BPData_nom.root"),opt.ptmin,opt.ptmax)
    sig,cuts,signeff,backeff=utils.significancecalc(data_x,100,variable,fs,fb)
    
    plt.figure()
    plt.clf()
    plt.plot(cuts,sig,label= f' max sig: {max(sig):.3f}; ' + f'cut:{cuts[np.argmax(sig)]:.2f};')
    plt.xlabel(variable +' cut')
    plt.ylabel('Significance')
    plt.legend()
    plt.savefig("../results/train/%s/significance.pdf" % config , bbox_inches='tight')

    plt.clf()
    plt.xlabel("DNN_output cut")
    plt.ylabel("Efficiency (Purity)")
    plt.plot(cuts, signeff,label="Signal")
    plt.plot(cuts, backeff,label="background")
    plt.legend()
    plt.savefig('../results/train/%s/sig_back_eff.pdf' % config, bbox_inches='tight')


    plt.clf()
    corr=signal[fullstagename].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (signal)")
    plt.xticks(range(len(fullstagename)),fullstagename,rotation=45)
    plt.yticks(range(len(fullstagename)),fullstagename)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig("../results/train/corrmatrix_sig_{}_{}.pdf".format(opt.ptmin,opt.ptmax), bbox_inches='tight')

    
    plt.clf()
    corr=background[fullstagename].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (background)")
    plt.xticks(range(len(fullstagename)),fullstagename,rotation=45)
    plt.yticks(range(len(fullstagename)),fullstagename)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig("../results/train/corrmatrix_back_{}_{}.pdf".format(opt.ptmin,opt.ptmax), bbox_inches='tight')

    data_imp=pd.DataFrame(test_X.numpy())
    data_imp.columns=stage
    tab=[]
    #tab2=[]
    #ynom=predict(model,test_X)
    for char in stage:
        pred=0
        for i in range(10):
            data_shuf=data_imp.copy()
            random.shuffle(data_shuf[char])
            pred+=torch.nn.functional.softmax(model(torch.tensor(data_shuf.values)),dim=1)[:,1]
        pred=pred/10
        tab.append(abs(yscore-pred).mean().item())

        #pred2=predict(model,torch.tensor(data_shuf.values))
        #tab2.append((pred2==ynom).sum().item()/len(ynom))

    plt.clf()
    plt.title("Feature importance")
    plt.bar(stage,tab)
    plt.ylabel("DNN_output difference")
    plt.xticks(rotation=45)
    plt.savefig("../results/train/%s/feature_importance.pdf" % config, bbox_inches='tight')

    #plt.clf()
    #plt.title("Feature importance")
    #plt.bar(stage,tab2)
    #plt.ylabel("y_pred difference")
    #plt.xticks(rotation=45)
    #plt.savefig("results/%s/feature_importance_alt.pdf" % config, bbox_inches='tight')

if __name__ == '__main__':
    main()
    print("Done!")
