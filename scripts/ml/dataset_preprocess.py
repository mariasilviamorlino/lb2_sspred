'''
A class for preprocessing data, both for GOR and SVM training, via Database objects
Data structure for the dataset: list of lists with the format [profile, sec_struct, id, gor_pred, svm_pred]

Methods:
* Database initialization (starting from IDs and profile/dssp folders)
* Svm/Libsvm input generation
* Addition of libsvm predictions in the correpsondent dataset entries

The Database objects have the following structure:
self.db --> list of lists
self.legenda --> a dict showing to which info each position in the inner lists corresponds
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
                        # otherwise, also consider the empty profile (consider points that map to the origin as points to be predicted)
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


    def svmpred(self, svmpred_path):
        '''Modifies the database object in-place
        1st argument: path to the file containing predictions'''
        mapping = {'1': 'H', '2': 'E', '3': '-'}
        svmpred_file = open(svmpred_path)
        for entry in self.db:
            sspred = ''
            i = 0
            while i < len(entry[1]):
                chpred = svmpred_file.readline().strip()
                sspred = sspred + mapping[chpred]
                i += 1
            entry[4] = sspred
        svmpred_file.close()


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
