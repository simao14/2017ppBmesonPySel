
import argparse
import uproot
from matplotlib import pyplot as plt
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
    opt = parser.parse_args()
    
    if opt.meson == 0:
        stages = [0,2,4,7,8,9,11,15]
        meson = "BP"
        binmin = 3
        binmax = 60
    else: 
        stages = [0,1,2,3,4,5,6,7,8,9,11,15]
        meson = "Bs"
        binmin = 3
        binmax = 50

    print("files")
    fileS = uproot.open(f"/lstore/cms/simao/sample/{meson}MC_{binmin}_{binmax}_small.root")
    fileB = uproot.open(f"/lstore/cms/simao/sample/{meson}Data_{binmin}_{binmax}_small.root")
    
    #fileR = uproot.open(f"/user/s/smcosta/data/{meson}Data_nom.root")
    fileR=uproot.open(f"../results/rootfiles/{meson}Data_nom_BDT.root") #use this only for pt 3--5
    
    signal,background = utils.prepdata(fileS,fileB,opt.ptmin,opt.ptmax,opt.meson)
    
    stage=utils.varset(stages)
    signal_var=signal[stage]
    background_var=background[stage]
    
    x=pd.concat([signal_var,background_var],axis=0,ignore_index=True)
    
    y=pd.concat([signal["tag"],background["tag"]],ignore_index=True)

    dset = xgb.DMatrix(x.to_numpy(), label=y.to_numpy())
    
    # initialize the model
    
    bst = xgb.Booster()
    bst.load_model(f"../results/{meson}/models/Xgb_{binmin}_{binmax}.model")

    ypred=bst.predict(dset)


    stagelist=utils.replacespecial(stages)
    config = "BDT-{}-{}-{}".format(opt.ptmin, opt.ptmax,stagelist)

    if not os.path.exists(f"../results/{meson}/sel/{config}"):
        os.makedirs(f"../results/{meson}/sel/{config}")
    
    utils.varplot(stage,signal_var,background_var,config,"sel",meson)


    nn_fpr, nn_tpr, _ = metrics.roc_curve(y.to_numpy(), ypred)
    auc=metrics.auc(nn_tpr,1-nn_fpr)
    plt.clf()
    plt.xlabel("signal efficiency")
    plt.ylabel("background efficiency")
    plt.plot(nn_tpr, 1-nn_fpr, label=f'auc:{auc:.3f};')
    plt.legend()
    plt.grid(visible=True)
    plt.savefig(f'../results/{meson}/sel/{config}/ROC-Curve.pdf', bbox_inches='tight')
    
    data_x=pd.DataFrame(x.to_numpy())
    data_x.columns=stage
    data_x["label"]=y.to_numpy()
    data_x["BDT_output"]=ypred
    data_x["scale_weight"]=np.ones(len(y.to_numpy()))
    variable="BDT_output"
    print("getting factors")
    fs,fb = utils_factors.get_factors(fileS,fileR,opt.ptmin,opt.ptmax,opt.meson)
    sig,cuts,signeff,backeff=utils.significancecalc(data_x,100,variable,fs,fb)

    plt.figure()
    plt.clf()
    plt.plot(cuts,sig,label= f' max sig: {max(sig):.3f}; ' + f'cut:{cuts[np.argmax(sig)]:.2f};')
    plt.xlabel(variable +' cut')
    plt.ylabel('Significance')
    plt.legend()
    plt.savefig(f"../results/{meson}/sel/{config}/significance.pdf", bbox_inches='tight')

    plt.clf()
    plt.xlabel("BDT_output cut")
    plt.ylabel("Efficiency (Purity)")
    plt.plot(cuts, signeff,label="Signal")
    plt.plot(cuts, backeff,label="background")
    plt.legend()
    plt.savefig(f'../results/{meson}/sel/{config}/sig_back_eff.pdf', bbox_inches='tight')

    plt.clf()
    corr=signal[stage].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (signal)")
    plt.xticks(range(len(stage)),stage,rotation=45)
    plt.yticks(range(len(stage)),stage)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig(f"../results/{meson}/sel/corrmatrix_sig_{opt.ptmin}_{opt.ptmax}.pdf", bbox_inches='tight')

   
    plt.clf()
    corr=background[stage].corr()
    plt.matshow(corr)
    plt.title("Correlation matrix (background)")
    plt.xticks(range(len(stage)),stage,rotation=45)
    plt.yticks(range(len(stage)),stage)
    for (i,j),z in np.ndenumerate(corr):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color="white")
    plt.colorbar()
    plt.savefig(f"../results/{meson}/sel/corrmatrix_back_{opt.ptmin}_{opt.ptmax}.pdf", bbox_inches='tight')

    feat=pd.DataFrame.from_dict([bst.get_score(importance_type='gain')])
    feat.columns=stage
    plt.clf()
    plt.title("Feature importance")
    feat.T.plot.barh(legend=False)
    plt.xlabel("BDT_output difference")
    plt.savefig(f"../results/{meson}/sel/{config}/feature_importance.pdf", bbox_inches='tight')

    

if __name__ == '__main__':
    main()
    print("Done!")