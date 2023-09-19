import ROOT
from ROOT import TFile
from ROOT import TH1F
from ROOT import RooRealVar
from ROOT import RooDataHist
from ROOT import RooArgSet
from ROOT import RooDataSet
from ROOT import RooExponential
from ROOT import RooGaussian
from ROOT import RooArgList
from ROOT import RooAddPdf
from ROOT import RooFit
from ROOT import RooFormulaVar
from ROOT import RooProduct
import numpy as np
import math

def prepdata(fileS,fileB,ptmin, ptmax):
    
    treeS = fileS["ntKp"]
    treeB = fileB["ntKp"]
    
    signal = treeS.arrays(library="pd")
    background = treeB.arrays(library="pd")
    
    cutS = ( (signal.Bgen==23333) & (signal.Bpt>ptmin) & (signal.Bpt<ptmax) )
    cutB = ((background.Bpt>ptmin) & (background.Bpt<ptmax) )
    
    signal_cut=signal[cutS]
    background_cut=background[cutB]

    return signal_cut,background_cut

def get_factors(fileS,fileB,ptmin, ptmax):

    bins = 100

    mmin = 5.0
    mmax = 6.0
    dmc, dd = prepdata(fileS,fileB,ptmin, ptmax)

    hmass = TH1F("hmass",";B-candidate mass (GeV)", bins, mmin, mmax)

    for index, event in enumerate(dd.Bmass):
        hmass.Fill(event)

    mass = RooRealVar("mass", "B-candidate mass", mmin, mmax, "GeV")
    args = RooArgList(mass)
    dh = RooDataHist("dh", "dh", args, hmass)

    Lambda = RooRealVar("Lambda", "lambda", -1.5, -5.0, 1.0)
    background = RooExponential("background", "background", mass, Lambda)


    scale = RooRealVar("scale","scale",1,0,2)
    mean = RooRealVar("mean", "mean", 0.5*(mmin+mmax), mmin, mmax)
    sigma1 = RooRealVar("sigma1", "sigma1", 0.1*(mmax-mmin),0.01,0.1)
    sigma2 = RooRealVar("sigma2", "sigma2", 0.1*(mmax-mmin),0.01,0.1)

    scaled_sigma1 = RooProduct("scaled_sigma1","scaled_sigma1",RooArgList(scale,sigma1))
    scaled_sigma2 = RooProduct("scaled_sigma2","scaled_sigma2",RooArgList(scale,sigma2))

    signal1 = RooGaussian("signal1", "signal1", mass, mean, scaled_sigma1)
    signal2 = RooGaussian("signal2", "signal2", mass, mean, scaled_sigma2)

    n_signal_initial = 0.8*dh.sumEntries()
    n_back_initial = 0.2*dh.sumEntries()
    n_signal = RooRealVar("n_signal","n_signal",n_signal_initial,-1,dh.sumEntries())
    n_back = RooRealVar("n_back","n_back",n_back_initial,-1,dh.sumEntries())

    sig1_frac = RooRealVar("sig1_frac", "sig1_frac", 0.2, 0., 1.)

    signal = RooAddPdf("signal", "signal", RooArgList(signal1,signal2), sig1_frac)

    model = RooAddPdf("model", "model", RooArgList(signal, background), RooArgList(n_signal, n_back))

    model.fitTo(dh)

    Mean = mean.getVal()
    Sigma = math.sqrt(sig1_frac.getVal() * sigma1.getVal()**2 + (1 - sig1_frac.getVal()) * sigma2.getVal()**2)

    sideband_edge1 = Mean - 3.5 * Sigma
    sideband_edge2 = Mean + 3.5 * Sigma
    
    mass_factor = RooRealVar("mass_factor", "B-candidate mass", mmin, mmax, "GeV")

    mass_factor.setRange("bkg range 1", mmin, sideband_edge1)
    mass_factor.setRange("bkg range 2", sideband_edge2, mmax)
    mass_factor.setRange("bkg range 3", sideband_edge1, sideband_edge2)

    argset = RooArgSet(mass_factor)

    Lambda_fit = RooRealVar("Lambda_fit", "lambda_fit", Lambda.getVal(), "GeV^-1")
    background_factor = RooExponential("background_factor", "background_factor", mass_factor, Lambda_fit)

    n_back_val=n_back.getVal()
    n_sig_val=n_signal.getVal()

    bkg_integral1 = background_factor.createIntegral(argset, RooFit.NormSet(argset), RooFit.Range("bkg range 1"))
    bkg_integral2 = background_factor.createIntegral(argset, RooFit.NormSet(argset), RooFit.Range("bkg range 2"))
    bkg_integral3 = background_factor.createIntegral(argset, RooFit.NormSet(argset), RooFit.Range("bkg range 3"))

    bkg_int1 = bkg_integral1.getVal()
    bkg_int2 = bkg_integral2.getVal()
    bkg_int3 = bkg_integral3.getVal()

    #Obtain the number of background events in each range
    bkg_events1 = bkg_int1 * n_back_val
    bkg_events2 = bkg_int2 * n_back_val
    bkg_events3 = bkg_int3 * n_back_val

    fb = bkg_events3 / (bkg_events1+bkg_events2)

    fs = n_sig_val/len(dmc)

    return fs,fb