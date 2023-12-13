
import argparse
import uproot
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as xgb
import pandas as pd
import numpy as np
import os.path
from sklearn import metrics

import sys
sys.path.insert(0, '../')
import utils
import utils_factors
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('meson', type=int)
    parser.add_argument('ptmin', type=float)
    parser.add_argument('ptmax', type=float)
    parser.add_argument('-iter', default=20 , type=int)
    parser.add_argument('-lamb', default=0.0001, type=float)
    parser.add_argument('-alpha', default= 0.13, type=float)
    parser.add_argument('-max_depth', default=9, type=int)
    parser.add_argument('-eta', type=float, default=0.32)
    parser.add_argument('-gamma', type=float, default= 0.0011)
    parser.add_argument('-rate_drop', type=float, default=2.082362173704635e-08)
    parser.add_argument('-skip_drop', type=float, default=0.002779566275170389)
    parser.add_argument('-booster',
                        choices=['gbtree', 'dart'], default='dart')
    parser.add_argument('-grow_policy',
                        choices= ["depthwise", "lossguide"], default='depthwise')
    parser.add_argument('-sample_type',
                        choices= ["uniform", "weighted"], default='weighted')
    parser.add_argument('-normalize_type',
                        choices=["tree", "forest"], default='tree')
    opt = parser.parse_args()
    
    print("files")
    if opt.meson == 0:
        stages = [0,2,4,7,8,9,11,15]
        fullstages = [0,2,4,7,8,9,11,12,14,15]
        fullstages2 = [0,2,4,7,8,9,11,12,14,15,16]
        meson = "BP"
        binmin = 3
        binmax = 60
    else: 
        stages = [0,1,2,3,4,5,6,7,8,9,11,15]
        fullstages = [0,1,2,3,4,5,6,7,8,9,11,12,14,15]
        fullstages2 = [0,1,2,3,4,5,6,7,8,9,11,12,14,15,16]
        meson = "Bs"
        binmin = 3
        binmax = 50


    fileS = uproot.open(f"/lstore/cms/simao/sample/{meson}MC_{binmin}_{binmax}_small.root")
    fileB = uproot.open(f"/lstore/cms/simao/sample/{meson}Data_{binmin}_{binmax}_small.root")
    
    
    signal,background = utils.prepdata(fileS,fileB,opt.ptmin,opt.ptmax,opt.meson)
    
    stage=utils.varset(stages)
    signal_var=signal[stage]
    background_var=background[stage]

    fullstagename=utils.varset(fullstages)
    
    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    train_x, test_x, train_y, test_y = train_test_split(x.to_numpy(), y.to_numpy(), test_size=0.2)
    train_x, valid_x, train_y, valid_y = train_test_split(train_x, train_y, test_size=0.2)
    
    dtrain = xgb.DMatrix(train_x, label=train_y)
    dvalid = xgb.DMatrix(valid_x, label=valid_y)
    dtest = xgb.DMatrix(test_x, label=test_y)

    # initialize the model

    param = {
        "verbosity": 1,
        "objective": "binary:logistic",
        "eval_metric": ["auc", "error"],
        "booster": opt.booster,
        "lambda": opt.lamb,
        "alpha": opt.alpha,
        "max_depth": opt.max_depth,
        "eta": opt.eta,
        "gamma": opt.gamma,
        "grow_policy": opt.grow_policy,
    }
    
    if param["booster"] == "dart":
        param["sample_type"] = opt.sample_type
        param["normalize_type"] = opt.normalize_type
        param["rate_drop"] = opt.rate_drop
        param["skip_drop"] = opt.skip_drop

    evals_result={}
    epochs=np.arange(0,opt.iter,1)
    bst = xgb.train(param, dtrain, opt.iter , evals=[(dtrain,"train"),(dvalid, "validation")],evals_result=evals_result)

    if not os.path.exists(f"../results/{meson}/models"):
        os.makedirs(f"../results/{meson}/models")
    
    
    bst.save_model(f"../results/{meson}/models/Xgb_{int(opt.ptmin)}_{int(opt.ptmax)}.model")

    stagelist=utils.replacespecial(stages)
    config = "BDT-{}-{}-{}".format(opt.ptmin, opt.ptmax,stagelist)

    if not os.path.exists(f"../results/{meson}/train/{config}"):
        os.makedirs(f"../results/{meson}/train/{config}")
    
    utils.varplot(stage,signal_var,background_var,config,"train",meson)

    plt.clf()
    plt.plot(epochs,evals_result["train"]["error"],label="train error")
    plt.plot(epochs,evals_result["validation"]["error"],label="validation error")
    plt.legend()
    plt.xlabel("Iteration")
    plt.ylabel("Error")
    plt.savefig(f"../results/{meson}/train/{config}/accuracy_per_iter.pdf", bbox_inches='tight')

    plt.clf()
    plt.plot(epochs,evals_result["train"]["auc"],label="train AUC")
    plt.plot(epochs,evals_result["validation"]["auc"],label="validation AUC")
    plt.legend()
    plt.xlabel("Iteration")
    plt.ylabel("AUC")
    plt.savefig(f"../results/{meson}/train/{config}/ROC_auc_per_iter.pdf", bbox_inches='tight')
    
    yscore = bst.predict(dtest)
    nn_fpr, nn_tpr, _ = metrics.roc_curve(test_y, yscore)
    auc=metrics.auc(nn_tpr,1-nn_fpr)

    plt.clf()
    plt.xlabel("signal efficiency")
    plt.ylabel("background efficiency")
    plt.plot(nn_tpr, 1-nn_fpr, label=f'auc:{auc:.3f};')
    plt.legend()
    plt.grid(visible=True)
    plt.savefig(f'../results/{meson}/train/{config}/ROC-Curve.pdf', bbox_inches='tight')
    

    data_x=pd.DataFrame(test_x)
    data_x.columns=stage
    data_x["label"]=test_y
    data_x["BDT_output"]=yscore
    data_x["scale_weight"]=np.ones(len(test_y))
    variable="BDT_output"
    fs,fb = utils_factors.get_factors(fileS,uproot.open(f"/user/s/smcosta/data/{meson}Data_nom.root"),opt.ptmin,opt.ptmax,opt.meson)
    sig,cuts,signeff,backeff=utils.significancecalc(data_x,100,variable,fs,fb)
    
    plt.figure()
    plt.clf()
    plt.plot(cuts,sig,label= f' max sig: {max(sig):.3f}; ' + f'cut:{cuts[np.argmax(sig)]:.2f};')
    plt.xlabel(variable +' cut')
    plt.ylabel('Significance')
    plt.legend()
    plt.savefig(f"../results/{meson}/train/{config}/significance.pdf", bbox_inches='tight')

    plt.clf()
    plt.xlabel("BDT_output cut")
    plt.ylabel("Efficiency (Purity)")
    plt.plot(cuts, signeff,label="Signal")
    plt.plot(cuts, backeff,label="background")
    plt.legend()
    plt.savefig(f'../results/{meson}/train/{config}/sig_back_eff.pdf', bbox_inches='tight')



    

    plt.clf()
    corr=signal[fullstagename].corr()
    fig,ax=plt.subplots(figsize=(8,8))
    # Add a color bar
    cbar = plt.colorbar(ax.matshow(corr, cmap='viridis'), ax=ax, fraction=0.046, pad=0.04)
    #cbar.set_label('Color Scale',fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    ax.set_xticks(range(len(fullstagename)),fullstagename,fontsize=14,rotation=45,ha='right')
    ax.set_yticks(range(len(fullstagename)),fullstagename,fontsize=14)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_title('Correlation matrix (signal)',fontsize=16)
    for (i, j), z in np.ndenumerate(corr):
        color = 'black' if z >= 0.6 else 'white'
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color=color, fontsize=14)

    fig.tight_layout()
    plt.savefig(f"../results/{meson}/train/corrmatrix_sig_{opt.ptmin}_{opt.ptmax}.pdf", bbox_inches='tight')

   
    plt.clf()
    corr=background[fullstagename].corr()
    fig,ax=plt.subplots(figsize=(8,8))
    # Add a color bar
    cbar = plt.colorbar(ax.matshow(corr, cmap='viridis'), ax=ax, fraction=0.046, pad=0.04)
    #cbar.set_label('Color Scale',fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    ax.set_xticks(range(len(fullstagename)),fullstagename,fontsize=13,rotation=45,ha='right')
    ax.set_yticks(range(len(fullstagename)),fullstagename,fontsize=13)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_title('Correlation matrix (back)',fontsize=16)
    for (i, j), z in np.ndenumerate(corr):
        color = 'black' if z >= 0.6 else 'white'
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color=color, fontsize=13)

    fig.tight_layout()

    plt.savefig(f"../results/{meson}/train/corrmatrix_back_{opt.ptmin}_{opt.ptmax}.pdf", bbox_inches='tight')

    plt.clf()
    fullstagename2=utils.varset(fullstages2)
    corr=pd.concat([signal,background],axis=0,ignore_index=True)[fullstagename2].corr()
    fig,ax=plt.subplots(figsize=(8,8))
    # Add a color bar
    cbar = plt.colorbar(ax.matshow(corr, cmap='viridis'), ax=ax, fraction=0.046, pad=0.04)
    #cbar.set_label('Color Scale',fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    ax.set_xticks(range(len(fullstagename2)),fullstagename2,fontsize=13,rotation=45,ha='right')
    ax.set_yticks(range(len(fullstagename2)),fullstagename2,fontsize=13)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_title('Correlation matrix (signal + back)',fontsize=16)
    for (i, j), z in np.ndenumerate(corr):
        color = 'black' if z >= 0.6 else 'white'
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color=color, fontsize=13)

    fig.tight_layout()
    plt.savefig(f"../results/{meson}/train/corrmatrix_all_{opt.ptmin}_{opt.ptmax}.pdf", bbox_inches='tight')

    feat=pd.DataFrame.from_dict([bst.get_score(importance_type='gain')])
    feat.columns=stage
    plt.clf()
    plt.title("Feature importance")
    feat.T.plot.bar(legend=False)
    plt.ylabel("BDT_output difference")
    plt.xticks(rotation=45)
    plt.savefig(f"../results/{meson}/train/{config}/feature_importance.pdf", bbox_inches='tight')

if __name__ == '__main__':
    main()
    print("Done!")
