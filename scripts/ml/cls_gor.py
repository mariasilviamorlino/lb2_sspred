'''
In this implementation, the GOR model interacts with database objects as defined in dataset_preprocess.py'''

import numpy as np

indexing = 'HE-'



class Gor:
    def __init__(self, w=17):
        self.w = w
        self.res_h = np.zeros((w,20))
        self.res_e = np.zeros((w,20))
        self.res_c = np.zeros((w,20))


    def __str__(self):
        return 'Rk,H\n{0}\nRk,E\n{1}\nRk,C\n{2}'.format(self.res_h, self.res_e, self.res_c)

    
    def train(self, datab):
        '''1st argument: database object
        '''
        size = 0
        padding = np.zeros((self.w//2,20))
        for i in datab:
            # preprocess input profile and dssp
            pssm = i[0]
            dssp = i[1]    
            pad_pr = np.vstack((padding,pssm,padding))
            pad_ss = ' '*(self.w//2) + dssp + ' '*(self.w//2)
            p, q, r = 0, self.w//2, self.w
            # increment counters
            while r <= len(pad_ss):
                if pad_ss[q] == 'H':
                    self.res_h += pad_pr[p:r]
                elif pad_ss[q] == 'E':
                    self.res_e += pad_pr[p:r]
                elif pad_ss[q] == '-':
                    self.res_c += pad_pr[p:r]
                p, q, r = p+1, q+1, r+1
            size += len(dssp)
        # normalize counts for dataset size (transform counts in frequencies)
        # print("\n\nCounts:\n{0}".format(self)) OK
        self.res_h = self.res_h/size
        self.res_e = self.res_e/size
        self.res_c = self.res_c/size
        # use the frequencies to compute information values
        tot_res_count = self.res_h + self.res_e + self.res_c
        self.res_h = self.res_h / (tot_res_count*(self.res_h.sum()/tot_res_count.sum()))
        self.res_h = np.log(self.res_h)
        self.res_e = self.res_e / (tot_res_count*(self.res_e.sum()/tot_res_count.sum()))
        self.res_e = np.log(self.res_e)
        self.res_c = self.res_c / (tot_res_count*(self.res_c.sum()/tot_res_count.sum()))
        self.res_c = np.log(self.res_c)


    def predict(self, datab, outfile=False):
        '''1st argument: Database object
        modify the Database object in-place (return new db value)
        pad the profile
        for every window: multiply profile to self.h self.e self.c and select the max'''
        padding = np.zeros((self.w//2,20))
        predss = '-HE'
        if outfile: out = open(outfile, 'w')
        for i in datab:
            # preprocess input profile and dssp
            pssm = i[0]
            pad_pr = np.vstack((padding,pssm,padding))
            p, q, r = 0, self.w//2, self.w
            pred = ''
            while r <= len(pad_pr):
                pred_h = np.sum(pad_pr[p:r]*self.res_h)
                pred_e = np.sum(pad_pr[p:r]*self.res_e)
                pred_c = np.sum(pad_pr[p:r]*self.res_c)
                l = [pred_c, pred_h, pred_e]  # since we are maxing on this list, if there are ties they are broken in favour of C
                # append to the prediction string the letter corresponding to the highest information
                pred = pred + predss[l.index(max(l))]
                p, q, r = p+1, q+1, r+1
            i[3] = pred
            if len(pred) != len(i[1]):
                print("Error; check the database entry {0}".format(i[2]))
                raise SystemExit
            if outfile: out.write('>{0}\n{1}\n\n'.format(i[2], i[3]))
        return datab



