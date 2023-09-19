
import argparse
import uproot
from matplotlib import pyplot as plt
import torch
import torch.nn as nn
import random
import pandas as pd
import numpy as np
import os.path
from sklearn import metrics
import sys
sys.path.insert(0, '../')
import utils
import utils_factors

def plot(epochs, plottable, xlabel="Epochs", ylabel='', name='',label=None):
    plt.clf()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(epochs, plottable, label=label)
    plt.legend()
    plt.grid(visible=True)
    plt.savefig('../results/sel/%s.pdf' % (name), bbox_inches='tight')
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ptmin', type=float)
    parser.add_argument('ptmax', type=float)
    parser.add_argument('-stages', type=int, nargs='+', default=[0,2,4,7,8,9,11,15]) #[0,2,4,7,8,9,11,12,15] #0,2,4,7,8,11
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
    
    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    dataset = utils.ClassificationDataset(x.values,y.values)
    
    # initialize the model
    
    model = torch.jit.load('../results/models/Py_3.0_60.0.pt')
    model.eval()

    ypred = torch.nn.functional.softmax(model(dataset.X),dim=1)[:,1]


    stagelist=utils.replacespecial(opt.stages)
    config = "NN-{}-{}-{}".format(opt.ptmin, opt.ptmax,stagelist)

    if not os.path.exists("../results/sel/%s" % config):
        os.makedirs("../results/sel/%s" % config)
    
    utils.varplot(stage,signal_var,background_var,config,"sel")


    nn_fpr, nn_tpr, _ = metrics.roc_curve(y.to_numpy(), ypred.detach().numpy())
    auc=metrics.auc(nn_tpr,1-nn_fpr)
    plot(nn_tpr,1-nn_fpr,ylabel="background efficiency", xlabel="signal efficiency",name='{}/ROC-curve'.format(config),label=f'auc:{auc:.3f};')
    
    data_x=pd.DataFrame(x.to_numpy())
    data_x.columns=stage
    data_x["label"]=y.to_numpy()
    data_x["NN_output"]=ypred
    data_x["scale_weight"]=np.ones(len(y.to_numpy()))
    variable="NN_output"
    print("getting factors")
    fs,fb = utils_factors.get_factors(fileS,uproot.open("/user/s/smcosta/data/BPData_nom.root"),opt.ptmin,opt.ptmax)
    sig,cuts,signeff,backeff=utils.significancecalc(data_x,100,variable,fs,fb)

    plt.figure()
    plt.clf()
    plt.plot(cuts,sig,label= f' max sig: {max(sig):.3f}; ' + f'cut:{cuts[np.argmax(sig)]:.2f};')
    plt.xlabel(variable +' cut')
    plt.ylabel('Significance')
    plt.legend()
    plt.savefig("../results/sel/%s/significance.pdf" % config , bbox_inches='tight')

    plt.clf()
    plt.xlabel("NN_output cut")
    plt.ylabel("Efficiency (Purity)")
    plt.plot(cuts, signeff,label="Signal")
    plt.plot(cuts, backeff,label="background")
    plt.legend()
    plt.savefig('../results/sel/%s/sig_back_eff.pdf' % config, bbox_inches='tight')

    plt.clf()
    corr=signal[stage].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (signal)")
    plt.xticks(range(len(stage)),stage,rotation=45)
    plt.yticks(range(len(stage)),stage)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig("../results/sel/corrmatrix_sig_{}_{}.pdf".format(opt.ptmin,opt.ptmax), bbox_inches='tight')

   
    plt.clf()
    corr=background[stage].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (background)")
    plt.xticks(range(len(stage)),stage,rotation=45)
    plt.yticks(range(len(stage)),stage)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig("../results/sel/corrmatrix_back_{}_{}.pdf".format(opt.ptmin,opt.ptmax), bbox_inches='tight')


    data_imp=pd.DataFrame(dataset.X.detach().numpy())
    data_imp.columns=stage
    tab=[]
    for char in stage:
        pred=0
        for i in range(10):
            data_shuf=data_imp.copy()
            random.shuffle(data_shuf[char])
            pred+=torch.nn.functional.softmax(model(torch.tensor(data_shuf.values)),dim=1)[:,1]
        pred=pred/10
        tab.append(abs(ypred-pred).mean().item())


    plt.clf()
    plt.title("Feature importance")
    plt.bar(stage,tab)
    plt.ylabel("DNN_output difference")
    plt.xticks(rotation=45)
    plt.savefig("../results/sel/%s/feature_importance.pdf" % config, bbox_inches='tight')

    

if __name__ == '__main__':
    main()
    print("Done!")