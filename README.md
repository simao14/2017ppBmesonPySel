# 2017ppBmesonPySel

## Setup

To be able to run this codes the first thing one should do is get their data ready for it. The codes in this git work based on python and therefore we will need to pass our data to numpy arrays. Normally root files are too big to pass all the events contained in them without numpy crashing. For this reason you can find the smallfile.C macro as an example on how to reduce the dimensions of your data. Try and only choose the variables you will be using on the analysis and do the most amount of cuts you can on this step for an easier time going forward. The macro here is only an example as it will be needed to adapt it for different analysis but one thing that should remain constant is that the final trees should be flat, that is an entry in a tree should be a single value, not an array.

Finally befora doing anything check the utils and utils_factor files. These have an assortment of functions useful throught the process and it would be best for the user to be familiar with them, especially since some of them apply cuts that might need changing in accord to the users needs.

## Optuna

We start in the train folder. First we find the optimal models utilizing optuna. OptBDTClassification.py is for the BDT models and OptClassification.py is for the NN models. You can modify the parameters of this code to better suit your needs. Note that you might need to change the varset function in the utils codes to correctly give the name of your variables. For brevety sake I will simply point to Optunas excelent documentation for further information. One thing I shall point out is that although the code gives the best performing model on the terminal I recommend installing and using optuna-dashboard to read the db file produced by this code as it gives a very in depth analysis of the optimization itself and can help guide a better choice of parameter ranges.

## training

Here you can set the default parameter values as the optimal ones obtained with optuna or call them when running the function. This macro will give you the saved trained model in /results/meson/models folder which will be used later for cutting the data. It also produces a variety of plots relative to the training which you can check out. Like the name implies XbgClassification is related to Xgboost (BDT) and PyClassification to Pytorch (NN).

## WP
Here we now go to the GetWorkingPoint folder.

This code, {Xgb/Py}GetPoint.py, is a simpler version of the previous one since its function is mainly to apply the saved models. Here you dont need to define any parameters since the models are already saved, just make sure the paths are right. All the functions related to the significance optimization are found in the utils_factor file. After the code is done you can check the significance plot and see the cut that maximizes it for the range given.

## Selection

Finally we apply the models and the cuts to the data itself, using the {Xgb/Py}Selection.py files. Here you have the usemc option. If zero you apply the selection to data, and if 1 it is applied to MC. You will need to do both for your analysis most likey. Change the cut values for what you obtained in the previous step. Like usual you might need to change the size of the array to fit your case but the ordering of the cuts should be in the same order as the ordering of the bin (if your bins are written smallest to highest then the cuts should be written from the one of the lowest bin to the highest obviously). After the models are used the code performs cuts in the data according to the WP values given and you obtain the final files with the selection applied.
