'''A set of functions to evaluate the performance of a prediction, given the predicted output and the expected one
Scoring indexes:
- (binary) MCC
- (binary) sensitivity
- (binary) positive predictive value
- (multi) three-class accuracy
- (segment-level, binary) SOV

SOV is defined as a separate function; all the confusion matrix-related indices are incorporated in the class confusion_matrix
'''

import numpy as np
import pandas as pd
import math
indexing = {'H': 0, 'E': 1, '-': 2}

class confusion_matrix:
    def __init__(self):
        self.cm = np.zeros((3, 3))
        self.bin_cm = None

    
    def populate(self, ypred, yreal):
        '''ypred and yreal are lists of the same length'''
        for i in range(len(ypred)):
            self.cm[indexing[ypred[i]], indexing[yreal[i]]]
    

    def binary_cm(self):
        '''return dictionary of all binary cms as numpy array, format: [0,0] TP; [0,1] FP; [1,0] FN; [1,1] TN'''
        bin_cm = {i:np.zeros((2,2)) for i in indexing.keys()}
        # populating H cm
        bin_cm['H'][0,0] = self.cm[0,0]
        bin_cm['H'][0,1] = self.cm[0,1:].sum()
        bin_cm['H'][1,0] = self.cm[1:,0].sum()
        bin_cm['H'][1,1] = self.cm[1:, 1:].sum()
        # populating E cm
        bin_cm['E'][0,0] = self.cm[1,1]
        bin_cm['E'][0,1] = self.cm[1,0] + self.cm[1,2]
        bin_cm['E'][1,0] = self.cm[0,1] + self.cm[2,1]
        bin_cm['E'][1,1] = self.cm[0,0] + self.cm[2,2] + self.cm[0,2] + self.cm[2,0]
        # populating C cm
        bin_cm['-'][0,0] = self.cm[2,2]
        bin_cm['-'][0,1] = self.cm[2,0:2].sum()
        bin_cm['-'][1,0] = self.cm[0:2,2].sum()
        bin_cm['-'][1,1] = self.cm[0:2,0:2].sum()
        self.bin_cm = bin_cm
        return bin_cm
    

    def sen(self, ss):
        '''1st argument: conformation with respect to which you want to compute the sensitivity'''
        bcm = self.bin_cm[ss]
        res = bcm[0,0]/(bcm[0,0]+bcm[1,0])
        return res
    

    def ppv(self, ss):
        '''1st argument: conformation with respect to which you want to compute the positive predictive value'''
        bcm = self.bin_cm[ss]
        res = bcm[0,0]/(bcm[0,0]+bcm[0,1])
        return res


    def mcc(self,ss):
        '''1st argument: conformation with respect to which you want to compute the mcc'''
        bcm = self.bin_cm[ss]
        tp, fp, fn, tn = bcm[0,0], bcm[0,1], bcm[1,0], bcm[1,1]
        res = (tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        return res
    

    def q3(self):
        res = (self.cm[0,0]+self.cm[1,1]+self.cm[2,2])/self.cm.sum()
        return res


    def performance_report(self):
        '''print a report formatted like a pandas dataframe'''
        pass


def sov(ypred, yreal):
    
    pass

# TODO:
# method for a performance report
