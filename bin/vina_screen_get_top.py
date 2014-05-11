#! /usr/bin/env python

import sys
import glob

def doit(n):
    file_names = glob.glob('*.pdbqt')
    everything = []
    for file_name in file_names:
        file = open(file_name)
        lines = file.readlines()
        file.close()
        line = lines[1]
        result = float(line.split(':')[1].split()[0])
        everything.append([result, file_name])
    everything.sort(lambda x,y: cmp(x[0], y[0]))
    part = everything[:n]
    for p in part:
        print p[1]

if __name__ == '__main__':
    doit(int(sys.argv[1]))
