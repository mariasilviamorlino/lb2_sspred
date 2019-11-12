'''paths
/home/mary/bioinfo/LB2/project/jpred4.list.txt
/home/mary/bioinfo/LB2/project/dssp
/home/mary/bioinfo/LB2/project/profiles/jpred_profiles
GOR method interacts with database objects as defined in dataset_preprocess.py'''

import numpy as np

indexing = 'HE-'


# def pad(obj, windowsize):
#     padding = np.zeros((windowsize//2,20))
#     return np.vstack((padding,obj,padding))


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
                size += len(dssp)
                p, q, r = p+1, q+1, r+1
        # normalize counts for dataset size
        self.res_h = self.res_h/size
        self.res_e = self.res_e/size
        self.res_c = self.res_c/size
        # transform the normalized counts in information values
        tot_res_count = self.res_h + self.res_e + self.res_c
        self.res_h = self.res_h / (tot_res_count*self.res_h.sum())
        self.res_h = np.log(self.res_h)
        self.res_e = self.res_e / (tot_res_count*self.res_e.sum())
        self.res_e = np.log(self.res_e)
        self.res_c = self.res_c / (tot_res_count*self.res_c.sum())
        self.res_c = np.log(self.res_c)
#        return self


    def predict(self, datab, outfile=False):
        '''1st argument: Database object
        modify the Database object in-place (return new db value)
        pad the profile
        for every window: multiply profile to self.h self.e self.c and select the max'''
        padding = np.zeros((self.w//2,20))
        predss = 'HE-'
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
            if outfile: out.write('>{0}\n{1}\n\n'.format(i[2], i[3]))
        return datab



if __name__ == '__main__':
    model = Gor()

    cavia = np.loadtxt('/home/mary/bioinfo/LB2/project/profiles/jpred_profiles/d3ag3g_.fasta.pssm')
    with open('/home/mary/bioinfo/LB2/project/dssp/d3ag3g_.dssp') as dsspfile:
        dsspseq = dsspfile.readlines()[1].strip()
#    dsspseq = 8*' '+dsspseq+8*' '
#    cavia = np.vstack((np.zeros((8,20)), cavia, np.zeros((8,20))))
#    print(model)

#     with open('/home/mary/bioinfo/LB2/project/jpred4.list.txt') as idlist:
#         for jid in idlist:
#             jid = jid.strip()
#             try:
#                 pr = np.loadtxt('/home/mary/bioinfo/LB2/project/profiles/jpred_profiles/'+jid+'.fasta.pssm')
#             except:
#                 pass
#             else:
# #                print(pr.shape)
#                 ss = open('/home/mary/bioinfo/LB2/project/dssp/'+jid+'.dssp').readlines()[1].strip()
#                 model.train(pr,ss)
# #                print(pr,ss)

    cavia2 = np.loadtxt('/home/mary/bioinfo/LB2/project/profiles/blind_profiles/4UFQ_A.fasta.pssm')

# TODO:
# fix the predict method so that it takes a Database object too
# test the whole thing and train the model in cross validation



