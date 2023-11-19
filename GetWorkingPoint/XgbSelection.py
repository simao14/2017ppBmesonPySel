
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('meson' , type=int)
    parser.add_argument('usemc' , type=int)
    opt = parser.parse_args()
    
    print("files")
    if opt.meson == 0:
        stages = [0,2,4,7,8,9,11,15]
        bin = [3,5,7,10,15,20,30,50,60]
        cut = [0.86,0.91,0.95,0.95,0.9,0.77,0.95,0.97]
        meson = "BP"
        binmin = 3
        binmax = 60
    else: 
        stages = [0,1,2,3,4,5,6,7,8,9,11,15]
        bin = [3,5,7,10,15,20,50]
        cut = [0.99,0.97,0.99,0.99,0.99,0.991]
        meson = "Bs"
        binmin = 3
        binmax = 50

    if opt.usemc == 0 :
        fileD = uproot.open(f"/lstore/cms/simao/sample/{meson}Data_{binmin}_{binmax}_small.root")
    else:
        fileD = uproot.open(f"/lstore/cms/simao/sample/{meson}MC_{binmin}_{binmax}_small.root")
    
    bst = xgb.Booster()
    bst.load_model(f"../results/{meson}/models/Xgb_{binmin}_{binmax}.model")
    
    stage=utils.varset(stages)
    
    data_pred=pd.DataFrame()
   
    for i in range(len(bin)-1):

        print(f"working on {bin[i]}-{bin[i+1]}")
        data= utils.seldata(fileD,bin[i],bin[i+1],opt.meson)
        data_var=data[stage]
        dset = xgb.DMatrix(data_var.to_numpy())
        ypred=bst.predict(dset)
        data["score"]=ypred
        data["BDT_pt_5_7"]=ypred
        data["BDT_pt_7_10"]=ypred
        data["BDT_pt_10_15"]=ypred
        data["BDT_pt_15_20"]=ypred
        data["BDT_pt_20_50"]=ypred
        data["BDT_pt_50_60"]=ypred
        data_pred=pd.concat([data_pred,data],axis=0,ignore_index=True)

    if opt.meson == 0:
        BDT_cut=(((data_pred.Bpt>3) & (data_pred.Bpt<5) & (data_pred.score>cut[0])) | ((data_pred.Bpt>5) & (data_pred.Bpt<7) & (data_pred.score>cut[1])) | ((data_pred.Bpt>7) & (data_pred.Bpt<10) & (data_pred.score>cut[2])) | ((data_pred.Bpt>10) & (data_pred.Bpt<15) & (data_pred.score>cut[3])) | ((data_pred.Bpt>15) & (data_pred.Bpt<20) & (data_pred.score>cut[4])) | ((data_pred.Bpt>20) & (data_pred.Bpt<30) & (data_pred.score>cut[5])) | ((data_pred.Bpt>30) & (data_pred.Bpt<50) & (data_pred.score>cut[6])) | ((data_pred.Bpt>50) & (data_pred.Bpt<60) & (data_pred.score>cut[7])))
    else:
        BDT_cut=(((data_pred.Bpt>3) & (data_pred.Bpt<5) & (data_pred.score>cut[0])) | ((data_pred.Bpt>5) & (data_pred.Bpt<7) & (data_pred.score>cut[1])) | ((data_pred.Bpt>7) & (data_pred.Bpt<10) & (data_pred.score>cut[2])) | ((data_pred.Bpt>10) & (data_pred.Bpt<15) & (data_pred.score>cut[3])) | ((data_pred.Bpt>15) & (data_pred.Bpt<20) & (data_pred.score>cut[4])) | ((data_pred.Bpt>20) & (data_pred.Bpt<50) & (data_pred.score>cut[5])) )

    
    data_final=data_pred[BDT_cut]
    
    if not os.path.exists(f"../results/rootfiles"):
        os.makedirs(f"../results/rootfiles")

    print("saving to root file")

    if opt.usemc == 0 :
        newfile=uproot.recreate(f"../results/rootfiles/{meson}Data_nom_BDT.root")
    else:
        newfile=uproot.recreate(f"../results/rootfiles/{meson}MC_nom_BDT.root")
    
    TreeD=data_final.to_dict(orient="list")

    if opt.meson == 0:
        newfile["ntKp"]=TreeD
    else: 
        newfile["ntphi"]=TreeD

if __name__ == '__main__':
    main()
    print("Done!")