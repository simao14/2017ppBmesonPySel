
import argparse
import uproot
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as xgb
import pandas as pd
import utils
import numpy as np
import os.path
from sklearn import metrics
import sys
sys.path.insert(0, '../')
import utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('usemc' , type=int)
    parser.add_argument('-bins', type=int, nargs='+', default=[3,5,7,10,15,20,30,50,60])
    parser.add_argument('-stages', type=int, nargs='+', default=[0,2,4,7,8,9,11,15]) 
    parser.add_argument('-cuts', type=int, nargs='+', default=[0.5,0.91,0.95,0.95,0.9,0.77,0.95,0.97]) 
    opt = parser.parse_args()
    
    print("files")
    if opt.usemc == 0 :
        fileD = uproot.open("/lstore/cms/simao/sample/BPData_3_60_small2.root")
    else:
        fileD = uproot.open("/lstore/cms/simao/sample/BPMC_3_60_small2.root")
    
    bst = xgb.Booster()
    bst.load_model("../results/models/Xgb_3.0_60.0.model")
    bin=opt.bins
    cut=opt.cuts
    stage=utils.varset(opt.stages)
    
    data_pred=pd.DataFrame()
   
    for i in range(len(bin)-1):

        print(f"working on {bin[i]}-{bin[i+1]}")
        data= utils.seldata(fileD,bin[i],bin[i+1])
        data_var=data[stage]
        dset = xgb.DMatrix(data_var.to_numpy())
        ypred=bst.predict(dset)
        data["score"]=ypred
        data_pred=pd.concat([data_pred,data],axis=0,ignore_index=True)

    BDT_cut=(((data_pred.Bpt>3) & (data_pred.Bpt<5) & (data_pred.score>cut[0])) | ((data_pred.Bpt>5) & (data_pred.Bpt<7) & (data_pred.score>cut[1])) | ((data_pred.Bpt>7) & (data_pred.Bpt<10) & (data_pred.score>cut[2])) | ((data_pred.Bpt>10) & (data_pred.Bpt<15) & (data_pred.score>cut[3])) | ((data_pred.Bpt>15) & (data_pred.Bpt<20) & (data_pred.score>cut[4])) | ((data_pred.Bpt>20) & (data_pred.Bpt<30) & (data_pred.score>cut[5])) | ((data_pred.Bpt>30) & (data_pred.Bpt<50) & (data_pred.score>cut[6])) | ((data_pred.Bpt>50) & (data_pred.Bpt<60) & (data_pred.score>cut[7])))
    data_final=data_pred[BDT_cut]
    
    if not os.path.exists("../results/rootfiles"):
        os.makedirs("../results/rootfiles")

    print("saving to root file")

    if opt.usemc == 0 :
        newfile=uproot.recreate("../results/rootfiles/BPData_nom_BDT.root")
    else:
        newfile=uproot.recreate("../results/rootfiles/BPMC_nom_BDT.root")
    
    TreeD=data_final.to_dict(orient="list")
    newfile["ntKp"]=TreeD


if __name__ == '__main__':
    main()
    print("Done!")