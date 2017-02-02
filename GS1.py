# -*- coding: utf-8 -*-
"""
Genome Sequencing

Created on Wed Feb  1 11:48:22 2017

@author: GY
"""


# Print data in certain format.
def prtout(text):
    for i in text:
        print(i)



# Generating k-mer collection in lexicographic order from a sequence.
# This is a simple illustration - what is in need is its inverse.
def composition(k, genome):
    col = []
    l = len(genome)
    for i in range(l-k+1):
        col.append(genome[i:(i+k)])
    #col.sort()
    return col



# Simplest trial: assemble sequence from well-ordered k-mer patterns.
def GenomePathOrd(patterns):
    seq = patterns[0]
    for i in patterns[1:]:
        seq += i[-1]
    return seq

    

# Mapping overlaping patterns:
def Overlap(patterns):
    n = len(patterns)
    mapdict = dict()
    for i in range(n-1):
        pat = patterns[i]
        if (pat[1:] in mapdict.keys()) and (pat[:-1] in mapdict.values()):
            continue
        if not pat[:-1] in mapdict.values():
            pre = pat[:-1]
            for j in range(i+1, n):
                if patterns[j][1:] == pre:
                    mapdict[patterns[j]] = pat
                    break
        if not pat[1:] in mapdict.keys():
            suf = pat[1:]
            for j in range(i+1, n):
                if patterns[j][:-1] == suf:
                    mapdict[pat] = patterns[j]
                    break
    return mapdict



# Plot overlap graph connecting former -> latter.
def PlotOverlap(patterns):
    mapdict = Overlap(patterns)
    for key, val in mapdict.items():
        print(key, '->', val)



# Mapping overlaping patterns using De Bruijn method.
def DeBruijn(k, genome):
    nodecol = composition(k-1, genome)
    n = len(nodecol)
    DB = dict()
    for i in range(n-1):
        if not nodecol[i] in DB.keys():
            DB[nodecol[i]] = [nodecol[i+1]]
        else:
            DB[nodecol[i]].append(nodecol[i+1])
    return DB



# DeBruijn graph from arbituary pattern set, rather than genome.
def DeBruijnP(k, patterns):
    DB = dict()
    for pattern in patterns:
        pre = pattern[:-1]
        suf = pattern[1:]
        if not pre in DB.keys():
            DB[pre] = [suf]
        else:
            DB[pre].append(suf)
    return DB



# Plotting De Bruijn graph.
def PlotDB(k, genome):
    if type(genome) is list:
        DB = DeBruijnP(k, genome)
    else:
        DB = DeBruijn(k, genome)
    for key in DB.keys():
        print(key, ' -> ', DB[key][0], sep = '', end = '')
        if len(DB[key]) > 1:
            for i in DB[key][1:]:
                print(',', i, sep = '', end = '')
        print('\n', sep = '', end = '')



# Output DeBruijn Graph to .txt file.
def WriteDB(k, genome):
    if type(genome) is list:
        DB = DeBruijnP(k, genome)
    else:
        DB = DeBruijn(k, genome)
    with open('task.txt', 'w') as f:
        for key in DB.keys():
            f.write(key + ' -> ' + DB[key][0])
            if len(DB[key]) > 1:
                for i in DB[key][1:]:
                    f.write(',' + i)
            f.write('\n')

            
            
  