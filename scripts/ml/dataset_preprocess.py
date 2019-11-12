'''
A toolbox for preprocessing data, both for GOR and SVM training
Data structure for the dataset: list of lists with the format (profile, sec_struct)
profile = numpy array; ss = string

what is in the toolbox?
* a function that, given in input a list of IDs and the location of the profiles and sec structs, returns a database
* a function to transform the database in a suitable way for SVM input
* a function to pickle and load databases

The Database objects have the following structure:
self.db --> list of lists
self.legenda --> a dict to tell you to which info each position in the inner lists corresponds
self.indexing --> a dict that links the indices of the outer lists to the ids of the inner lists
'''

import numpy as np

# strings to keep indexing in mind
sec_struct = 'HE-'
amino_acids = 'ARNDCQEGHILKMFPSTWYV'
# mapping of the classes for multiclass svm
class_mapping = {'H': 1, 'E': 2, '-': 3}

class Database:
    '''
    Attributes:
    db --> list of lists [numpy_profile, string_ss, seq_id, ss_pred]'''

    def __init__(self, ids, profile_path, ss_path, log = False):
        '''1st argument: list of ids (list or file handler)
        2nd argument: path to pssm files
        3rd argument: path to dssp files
        kwarg log (boolean): generate output file with a list of all ids for which there is no pssm (default False)'''
        self.db = list()
        self.indexing = dict()
        self.legenda = {0: 'profile',
                        1: 'ss_real',
                        2: 'id',
                        3: 'ss_pred_GOR',
                        4: 'ss_pred_SVM'}
        ctr = 0
        if log: logfile = open('unavailable_profiles.txt', 'a')
        for i in ids:
            j = i.strip()
            try:
                prof = np.loadtxt(profile_path+'/'+j+'.fasta.pssm')
            except:
                if log:
                    logfile.write('{0}.fasta.pssm in {1}: not found.\n'.format(j, profile_path))
            else:
                prof = prof/100
                ssfile = open(ss_path+'/'+j+'.dssp')
                ss = ssfile.readlines()[1].strip()
                self.db.append([prof, ss, j, None, None])
                self.indexing[j] = ctr
                ssfile.close()
                ctr += 1
        if log: logfile.close()
    

    def __iter__(self):
        '''Note: the iteration only goes over self.db and neglects other attributes'''
        return dbIterator(self.db)


    def makesvminput(self, wsize = 17):
        '''
        Format the data in a suitable format for scikitlearn SVMs
        kwarg wsize (int): window size (necessary for padding and profile window vectorialization)
        '''
        # shape of the input:
        # examples are an array of 'vectorialized' profile windows
        # target classes are an array of characters referring to the central residue of the profile windows

        
        # output initialization
        x = list()
        y = list()
        for i in self.db:
            # padding
            padding = np.zeros((wsize//2,20))
            pad_pr = np.vstack((padding,i[0],padding))
            pad_ss = ' '*(wsize//2) + i[1] + ' '*(wsize//2)
            # window serialization
            p, q, r = 0, wsize//2, wsize
            while r <= len(pad_ss):
                xit = pad_pr[p]
                for n in range(p+1,r):
                    xit = np.hstack((xit,pad_pr[n]))
                x.append(xit)
                y.append(class_mapping[pad_ss[q]])
                p, q, r = p+1, q+1, r+1
        x = np.array(x)
        y = np.array(y)
        return x, y

    
    def makelibsvminput(self, filepath = 'training.dat.txt', wsize = 17, testing = False):
        '''
        kwarg filename: path of output file
        kwarg wsize: size of profile window to serialize
        kwarg testing: specifies whether the input is intended for training or for predicting (default False)
        return: none
        '''
        with open(filepath, 'w') as outfile:
            for i in self.db:
                # padding
                padding = np.zeros((wsize//2,20))
                pad_pr = np.vstack((padding,i[0],padding))
                pad_ss = ' '*(wsize//2) + i[1] + ' '*(wsize//2)

                # window serialization
                p, q, r = 0, wsize//2, wsize
                while r <= len(pad_ss):
                    xit = pad_pr[p]
                    for n in range(p+1,r):
                        xit = np.hstack((xit,pad_pr[n]))

                    if xit.sum() == 0 and not testing:
                        pass
                        # if the database is meant for training, write to file only if profile is not empty
                    else:
                        # otherwise, also consider the empty profile (predict a class for points that map to the origin)
                        line = str(class_mapping[pad_ss[q]])
                        ctr = 1
                        # this for cycle takes care of the sparse notation
                        for elem in xit:
                            if elem != 0:
                                line = '{0} {1}:{2}'.format(line, ctr, elem)
                            ctr += 1
                        outfile.write(line+'\n')
                    # increment counters
                    p += 1; q += 1; r += 1
    

    def onehot(self, fastapath):
        '''Transforms zero matrices in one-hot encoded sequences. To be used when datasets are created for the purpose of testing'''
        for item in self.db:
            if item[0].sum() == 0:
                fastafile = open(fastapath+'/'+item[2]+'.fasta')
                fasta = fastafile.readlines()[1].strip()
                fastafile.close()
                ctr = 0
                for char in fasta:
                    # ctr indexes the position 
                    item[0][ctr, amino_acids.index(char)] = 1
                    ctr += 1


    def add_pred(self, seqid, ss_pred, SVM = True):
        '''
        adds a prediction made with GOR
        seqid: string - id of the sequence/profile/db entry
        ss_pred: string - prediction of secondary structure
        kwarg SVM: boolean, default true. If it's true, the prediction is added in the 4th slot of the list
        '''
        indx = self.indexing[seqid]
        self.db[indx][3] = ss_pred
        pass



class dbIterator:
    """Iterator class for Database objects; note: the iteration only goes over the list of lists!"""
    def __init__(self, datab):
        self.datab = datab
        self.index = 0
    
    def __next__(self):
        if self.index < len(self.datab):
            result = self.datab[self.index]
            self.index += 1
            return result
        raise StopIteration



if __name__ == '__main__':
    ids_path = '/home/mary/bioinfo/LB2/project/cv/fold1/cv1234.id'
    profiles_path = '/home/mary/bioinfo/LB2/project/profiles/jpred_profiles'
    dssp_path = '/home/mary/bioinfo/LB2/project/dssp'
    model_path = '/home/mary/bioinfo/LB2/project/models'

    idss = open(ids_path)

    testdb = Database(idss, profiles_path, dssp_path)
#    testdb.makelibsvminput(filepath = '/home/mary/bioinfo/LB2/project/libsvm/prova.dat.txt')

# TODO
# give the dataset a prediction attribute - same n° of items as self.db
