'''A set of functions to evaluate the performance of a prediction, given the predicted output and the expected one
Scoring indices:
- (binary) MCC
- (binary) sensitivity
- (binary) positive predictive value
- (multi) three-class accuracy

all the confusion matrix-related indices are implemented as methods of the class confusion_matrix
'''

import numpy as np
import math
indexing = {'H': 0, 'E': 1, '-': 2}

class confusion_matrix:
    def __init__(self):
        self.cm = np.zeros((3, 3))
        self.bin_cm = {i:np.zeros((2,2)) for i in indexing.keys()}

    
    def __str__(self):
        return "{0}\n\n{1}".format(self.cm, self.bin_cm)


    def populate(self, ypred, yreal):
        '''ypred and yreal are lists of the same length'''
        if len(ypred) != len(yreal):
            print(ypred, '\n\n\n', yreal)
            print("something is horribly wrong")
            raise SystemExit
        for i in range(len(ypred)):
            self.cm[indexing[ypred[i]], indexing[yreal[i]]] += 1
        self.binary_cm()
    

    def binary_cm(self):
        '''return dictionary of all binary cms as numpy array, format: [0,0] TP; [0,1] FP; [1,0] FN; [1,1] TN'''
        # populating H cm
        self.bin_cm['H'][0,0] = self.cm[0,0]
        self.bin_cm['H'][0,1] = self.cm[0,1:].sum()
        self.bin_cm['H'][1,0] = self.cm[1:,0].sum()
        self.bin_cm['H'][1,1] = self.cm[1:, 1:].sum()
        # populating E cm
        self.bin_cm['E'][0,0] = self.cm[1,1]
        self.bin_cm['E'][0,1] = self.cm[1,0] + self.cm[1,2]
        self.bin_cm['E'][1,0] = self.cm[0,1] + self.cm[2,1]
        self.bin_cm['E'][1,1] = self.cm[0,0] + self.cm[2,2] + self.cm[0,2] + self.cm[2,0]
        # populating C cm
        self.bin_cm['-'][0,0] = self.cm[2,2]
        self.bin_cm['-'][0,1] = self.cm[2,0:2].sum()
        self.bin_cm['-'][1,0] = self.cm[0:2,2].sum()
        self.bin_cm['-'][1,1] = self.cm[0:2,0:2].sum()
        return self.bin_cm
    

    def sen(self, ss):
        '''1st argument: conformation with respect to which you want to compute the sensitivity'''
        bcm = self.bin_cm[ss]
        res = bcm[0,0]/(bcm[0,0]+bcm[1,0])
        return 100*res
    

    def ppv(self, ss):
        '''1st argument: conformation with respect to which you want to compute the positive predictive value'''
        bcm = self.bin_cm[ss]
        res = bcm[0,0]/(bcm[0,0]+bcm[0,1])
        return 100*res


    def mcc(self,ss):
        '''1st argument: conformation with respect to which you want to compute the mcc'''
        bcm = self.bin_cm[ss]
        tp, fp, fn, tn = bcm[0,0], bcm[0,1], bcm[1,0], bcm[1,1]
        res = (tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        return 100*res
    

    def q3(self):
        cm = self.cm
        res = (cm[0,0]+cm[1,1]+cm[2,2])/cm.sum()
        return 100*res


if __name__ == '__main__':
    real = '-------HHHHHHHHHH---------EEEEEEE----EEEEEEE------'
    pred = '---HHHHHH-HHHHHHHHHHH--EEEEEEEEEEEE-------HHHEHE--'

    cmatrix = confusion_matrix()
    cmatrix.populate(pred,real)
    print(cmatrix)
    print(cmatrix.sen('H'))
    print(cmatrix.mcc('E'))
    print(cmatrix.q3())
