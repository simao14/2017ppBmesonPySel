import torch
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import os.path

class ClassificationDataset(torch.utils.data.Dataset):

    def __init__(self, data, labels):

        """
        data: the dict returned by utils.load_classification_data
        """
        
        train_X = data
        train_y = labels
        
        self.X = torch.tensor(train_X, dtype=torch.float32)
        self.y = torch.tensor(train_y, dtype=torch.long)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]
    
def replacespecial(stages):

    string=str(stages[0])
    for i in range(1,len(stages)):
        string+= "_"+ str(stages[i]) 
    return string

def varset(stages):
    
    var=["Btrk1Pt","Btrk2Pt","Trk1DCAz","Trk2DCAz","Trk1DCAxy","Trk2DCAxy","MassDis","dls","Balpha","dls2D","cos(Bdtheta)","Bchi2cl","Btrk1Eta","Btrk2Eta","Bmass","Bpt","tag"]
    stage=[]
    for i in range(len(stages)):
        stage.append(var[stages[i]])
    return stage 

def prepdata(fileS,fileB,ptmin, ptmax,meson_n):
    
    if meson_n==0:
      treeS = fileS["ntKp"]
      treeB = fileB["ntKp"]
      mass = 5.27929
      lside = 0.25
      rside = 0.30
    else:
      treeS = fileS["ntphi"]
      treeB = fileB["ntphi"]
      mass = 5.36682
      lside = 0.20
      rside = 0.30
    
    signal = treeS.arrays(library="pd")
    background = treeB.arrays(library="pd")

    signal["Trk1DCAz"] = abs(signal.Btrk1Dz1/signal.Btrk1DzError1)
    signal["Trk2DCAz"] = abs(signal.Btrk2Dz1/signal.Btrk2DzError1)
    signal["Trk1DCAxy"] = abs(signal.Btrk1Dxy1/signal.Btrk1DxyError1)
    signal["Trk2DCAxy"] = abs(signal.Btrk2Dxy1/signal.Btrk2DxyError1)
    signal["MassDis"] = abs(signal.Btktkmass-1.019455)
    signal["dls"] = signal.BsvpvDistance/signal.BsvpvDisErr
    signal["dls2D"] = signal.Bd0
    signal["tag"] = np.ones(signal.shape[0])
    
    background["Trk1DCAz"] = abs(background.Btrk1Dz1/background.Btrk1DzError1)
    background["Trk2DCAz"] = abs(background.Btrk2Dz1/background.Btrk2DzError1)
    background["Trk1DCAxy"] = abs(background.Btrk1Dxy1/background.Btrk1DxyError1)
    background["Trk2DCAxy"] = abs(background.Btrk2Dxy1/background.Btrk2DxyError1)
    background["MassDis"] = abs(background.Btktkmass-1.019455)
    background["dls"] = background.BsvpvDistance/background.BsvpvDisErr
    background["dls2D"] = background.Bd0
    background["tag"] = np.zeros(background.shape[0])
    
    
    cutS = ( (signal.Bgen==23333) & (signal.Bpt>ptmin) & (signal.Bpt<ptmax) )
    cutB = ( (((background.Bmass - mass ) > lside) &  ((background.Bmass - mass) < rside)) & (background.Bpt>ptmin) & (background.Bpt<ptmax) )
    
    signal_cut=signal[cutS]
    background_cut=background[cutB]

    return signal_cut,background_cut

def seldata(fileD,ptmin, ptmax,meson_n):
    
    if meson_n == 0:
      treeD = fileD["ntKp"]
    else:
       treeD = fileD["ntphi"]

    data = treeD.arrays(library="pd")
    
    data["Trk1DCAz"] = abs(data.Btrk1Dz1/data.Btrk1DzError1)
    data["Trk2DCAz"] = abs(data.Btrk2Dz1/data.Btrk2DzError1)
    data["Trk1DCAxy"] = abs(data.Btrk1Dxy1/data.Btrk1DxyError1)
    data["Trk2DCAxy"] = abs(data.Btrk2Dxy1/data.Btrk2DxyError1)
    data["MassDis"] = abs(data.Btktkmass-1.019455)
    data["dls"] = data.BsvpvDistance/data.BsvpvDisErr
    data["dls2D"] = data.Bd0
    
    cutS = ( (data.Bpt>ptmin) & (data.Bpt<ptmax) )
    
    data_cut=data[cutS]

    return data_cut

def significance (s,b):
  
  if math.sqrt(s+b)==0:
    sig=0
  else:
    sig = s/math.sqrt(s+b)
  return sig if not np.isnan(sig) else 0
 
def varplot(stage,signal,backround,outname,typename,meson):

  if not os.path.exists(f"../results/{meson}/{typename}/{outname}/Variables"):
        os.makedirs(f"../results/{meson}/{typename}/{outname}/Variables")

  for i in stage:
    plt.clf()
    plt.hist(signal[i],density=True,bins=100,label=f'Signal - nº entries:{len(signal[i])}; mean:{np.mean(signal[i]):.2f}; dev:{np.std(signal[i]):.2f};',color="r")
    plt.hist(backround[i],density=True,bins=100,label=f'Background - nº entries:{len(backround[i])}; mean:{np.mean(backround[i]):.2f}; dev:{np.std(backround[i]):.2f};',color="b" ,fill=False)
    plt.legend()
    plt.ylabel("nº of entries")
    plt.xlabel(i)
    plt.savefig(f"../results/{meson}/{typename}/{outname}/Variables/{i}_distribution.png")

def significancecalc(data,Ncuts,variable,fs,fb,quant=1):
  if quant==1:
    cuts = np.linspace(min(data[variable]),max(data[variable]),Ncuts) 
  else:
    cuts = np.linspace(min(data[variable]),data.query(f'{variable}< {variable}.quantile({quant})')[variable].max(),Ncuts) 
  sig  = [0 for i in range(Ncuts)]
  signeff  = [0 for i in range(Ncuts)]
  backeff  = [0 for i in range(Ncuts)]
  nsig=data.loc[data.query('(label==1)').index,'scale_weight'].sum()
  nback=data.loc[data.query('(label==0)').index,'scale_weight'].sum()
  for i,cut in enumerate(cuts):

    s = data.loc[data.query(f'(label==1) & ({variable}>{cut})').index,'scale_weight'].sum()
    b = data.loc[data.query(f'(label==0) & ({variable}>{cut})').index,'scale_weight'].sum()

    sig[i] = significance(s*fs,b*fb)
    signeff[i]=s/nsig
    backeff[i]=b/nback
    
  return sig, cuts, signeff, backeff