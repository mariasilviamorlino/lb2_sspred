def find_segments(ypred, yreal, ss):
    '''
    1st argument --> string of predicted secondary structure
    2nd argument --> string of observed secondary structure
    3rd argument --> helix, strand or coil
    ypred and yreal must have the same length'''
    # variable and flag initialization
    s1, s2, S1, S2 = list(), list(), list(), list()
    s1_ctr, s2_ctr = 0, 0
    map = list()
    s_h1 = list()
    flags1, flags2, mapped = False, False, False
    no_ov = True
    yprpad = ypred+'*'
    yrepad = yreal+'*'

    for i in range(len(yprpad)):
        # definition of s1 segments by start and end index
        if yrepad[i] == ss and yrepad[i] != '*':
            if not flags1:
                s1 = list()
                s1.append(i)
                flags1 = True
        else:
            if flags1:
                s1.append(i)
                flags1 = False
                S1.append(s1)
                # if the segment ends and no overlapping segment has been found,
                # add the segment to S'(H)
                if no_ov: s_h1.append(s1)
                else: no_ov = True  # reset the flag signalling the presence of overlaps
                s1_ctr += 1
        # definition of s2 segments by start and end index
        if yprpad[i] == ss and yprpad != '*':
            if not flags2:
                s2 = list()
                s2.append(i)
                flags2 = True
        else:
            if flags2:
                s2.append(i)
                flags2 = False
                S2.append(s2)
                s2_ctr += 1
        
        if yrepad[i] == ss and yprpad[i] == ss and yrepad[i] != '*':
            if not mapped:
                map.append((s1_ctr, s2_ctr))
                mapped = True
                no_ov = False
        else:
            # reset flags
            mapped = False

    return(S1, S2, map, s_h1)
    

def normalizer(SH, SH1):
    '''
    1st argument --> couples of overlapping segments. List of tuples of lists
    2nd argument --> non-overlapping yreal segments. List of lists'''
    norm = 0
    for i in SH:
        lens1 = i[0][1]-i[0][0]
        norm = norm + lens1
    for i in SH1:
        lens1 = i[1]-i[0]
        norm = norm + lens1
    return norm


def minov(s1, s2):
    start = [s1[0], s2[0]]
    end = [s1[1], s2[1]]
    minoverlap = min(end) - max(start)
    return minoverlap


def maxov(s1, s2):
    start = [s1[0], s2[0]]
    end = [s1[1], s2[1]]
    maxoverlap = max(end) - min(start)
    return maxoverlap


def delta(s1, s2):
    d = min(
        maxov(s1, s2) - minov(s1, s2),
        minov(s1, s2),
        (s1[1]-s1[0])//2,
        (s2[1]-s2[0])//2
    )
    return d


def summation(Sh):
    '''1st argument --> list of overlapping segments. List of tuples of lists
    S(H) = [ ([s1start, s1end], [s2start, s2end]) , ([..,..],[..,..]), (...)]'''
    sum = 0
    for couple in Sh:
        s1 = couple[0]
        s2 = couple[1]
        sum = sum + ((minov(s1,s2) + delta(s1,s2))/maxov(s1,s2))*(s1[1]-s1[0])
    return sum


def sov(ypr, yre, ss):
    '''1st argument --> predicted ss
    2nd argument --> real ss
    ss --> helix, strand or coil'''
    # scan the sequences to define S(H) and S'(H)
    # ypred and yreal are supposed to have the same length.
    # sh and sh1 are lists of tuples. They are the basis of SOV calculation
    S1, S2, map, sh1 = find_segments(ypr, yre, ss)
    sh = list()
    for couple in map:
        sh.append((S1[couple[0]], S2[couple[1]]))
    # return(sh, sh1)
    norm = normalizer(sh, sh1)
    sov_score = (summation(sh)*100)/norm
    return sov_score


def multiclass_sov(ypr, yre):
    norm = 0
    summ = 0
    for char in 'HE-':
        S1, S2, map, si1 = find_segments(ypr, yre, char)
        si = list()
        for couple in map:
            si.append((S1[couple[0]], S2[couple[1]]))
        norm = norm + normalizer(si, si1)
        summ = summ + summation(si)
    sov_tot = (summ * 100)/norm
    return sov_tot


if __name__ == '__main__':
    # testcode: sequences taken from the SOV'99 and SOV_refined articles, for which the reference score is known.
    ref = '-HHHHHHHHHH-'
    for i in (['-H-H-H-H-H--', 12.5], ['-HHH-HHH-HH-', 40.6], ['-HH--HHHHH--', 52.3], ['---HHHH-----', 54.4],
                ['---HHHHH----', 63.2], ['---HHHHHH---', 80.6], ['---HHHHHHH--', 90.3], ['---HHHHHHHH-', 94.4]):
        pred = i[0]
        print('my segment overlap: {0}\nreference value: {1}\n'.format(multiclass_sov(pred,ref), i[1]))
