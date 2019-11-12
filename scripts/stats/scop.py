import sys
from statstools import compo 

class_des = dict()
with open('/home/mary/bioinfo/LB2/project/scopclass.des.txt') as scopdes:
    for line in scopdes:
        l = line.split('\t')
        class_des[l[0].strip()] = l[1].strip()

if __name__ == '__main__':
    countrs = dict()
    clids = open('/home/mary/bioinfo/LB2/project/scopclass.txt')
    cllist = []
    for i in clids:
        cllist.append(i.strip())
    clids.close()
    stat = compo(cllist, countrs)
    print(stat)
    print(class_des)
