
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
    parser.add_argument('-iter', default=20 , type=int)
    parser.add_argument('-lamb', default=0.00010693059802121613, type=float)
    parser.add_argument('-alpha', default= 0.1272597633503785, type=float)
    parser.add_argument('-max_depth', default=9, type=int)
    parser.add_argument('-eta', type=float, default=0.32297132116252375)
    parser.add_argument('-gamma', type=float, default= 0.0011261590476240695)
    parser.add_argument('-rate_drop', type=float, default=2.082362173704635e-08)
    parser.add_argument('-skip_drop', type=float, default=0.002779566275170389)
    parser.add_argument('-stages', type=int, nargs='+', default=[0,2,4,7,8,9,11,15]) #[0,2,4,7,8,9,11,12,15] #0,2,4,7,8,11
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

    if not os.path.exists("../results/models"):
        os.makedirs("../results/models")
    
    
    bst.save_model(f"../results/models/Xgb_{opt.ptmin}_{opt.ptmax}.model")

    stagelist=utils.replacespecial(opt.stages)
    config = "BDT-{}-{}-{}".format(opt.ptmin, opt.ptmax,stagelist)

    if not os.path.exists("../results/train/%s" % config):
        os.makedirs("../results/train/%s" % config)
    
    utils.varplot(stage,signal_var,background_var,config,"train")

    plt.clf()
    plt.plot(epochs,evals_result["train"]["error"],label="train error")
    plt.plot(epochs,evals_result["validation"]["error"],label="validation error")
    plt.legend()
    plt.xlabel("Iteration")
    plt.ylabel("Error")
    plt.savefig(f"../results/train/{config}/accuracy_per_iter.pdf", bbox_inches='tight')

    plt.clf()
    plt.plot(epochs,evals_result["train"]["auc"],label="train AUC")
    plt.plot(epochs,evals_result["validation"]["auc"],label="validation AUC")
    plt.legend()
    plt.xlabel("Iteration")
    plt.ylabel("AUC")
    plt.savefig(f"../results/train/{config}/ROC_auc_per_iter.pdf", bbox_inches='tight')
    
    yscore = bst.predict(dtest)
    nn_fpr, nn_tpr, _ = metrics.roc_curve(test_y, yscore)
    auc=metrics.auc(nn_tpr,1-nn_fpr)
    plot(nn_tpr,1-nn_fpr,ylabel="background efficiency", xlabel="signal efficiency",name='{}/ROC-curve'.format(config),label=f'auc:{auc:.3f};')

    data_x=pd.DataFrame(test_x)
    data_x.columns=stage
    data_x["label"]=test_y
    data_x["BDT_output"]=yscore
    data_x["scale_weight"]=np.ones(len(test_y))
    variable="BDT_output"
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
    plt.xlabel("BDT_output cut")
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

    plt.clf()
    fullstagename2=utils.varset([0,2,4,7,8,9,11,12,14,15,16])
    corr=pd.concat([signal,background],axis=0,ignore_index=True)[fullstagename2].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (background + signal)")
    plt.xticks(range(len(fullstagename2)),fullstagename2,rotation=45)
    plt.yticks(range(len(fullstagename2)),fullstagename2)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig("../results/train/corrmatrix_all_{}_{}.pdf".format(opt.ptmin,opt.ptmax), bbox_inches='tight')

    feat=pd.DataFrame.from_dict([bst.get_score(importance_type='gain')])
    feat.columns=stage
    plt.clf()
    plt.title("Feature importance")
    feat.T.plot.barh(legend=False)
    plt.xlabel("BDT_output difference")
    plt.savefig("../results/train/%s/feature_importance.pdf" % config, bbox_inches='tight')

if __name__ == '__main__':
    main()
    print("Done!")
