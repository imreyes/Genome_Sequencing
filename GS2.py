# -*- coding: utf-8 -*-
"""
Genome Sequencing 2

Created on Thu Feb  2 14:03:40 2017

@author: GY
"""

import GS1
import copy as cp


# Read De Bruijn graph into dictionary.
def ReadDB(graph):
    if type(graph) is dict:             # What we want.
        return graph
    if type(graph) is str:              # First parse string to list.
        graph = graph.splitlines()
    if type(graph) is list:             # Then process list to dict.
        DB = dict()
        for item in graph:
            nums = item.split()
            DB[nums[0]] = nums[2].split(',')
    return DB


   
 
# Form Eulerian cycle from DB graph.
def EulerianCycle(DB):
    DB = ReadDB(DB)                     # Process DB into dict; allows 'raw' input.
    init = list(DB.keys())[0]
    totalcycle = []
    idx = 0
    while DB:                           # Cycle until DB empties - all paths covered.
        cycle = []
        while DB[init]:                 # Small cycle ends when there's no place to continue (no values with the given key).
            new = DB[init].pop()
            cycle.append(new)
            init = new
        DB = dict((k, v) for k, v in DB.items() if v)       # Clear off empty items.
        if totalcycle:                  # Concatenate subcycle into entire cycle.
            totalcycle = totalcycle[:idx] + cycle + totalcycle[idx:]
        else:
            totalcycle = cycle
        if DB:                          # Select one point in both remaining DB and the cycle.
            for i in DB.keys():
                if i in totalcycle:
                    init = i
                    idx = totalcycle.index(init) + 1
                    break
    return totalcycle



# Form Eulerian cycle from DB graph.
def EulerianPath(DB):
    DB = ReadDB(DB)
    DBCounter = dict()      # This dict stores numbers of in-paths of items.
    missedkey = None
    for key, val in DB.items():
        if not key in DBCounter.keys():
            DBCounter[key] = 0
        for j in val:
            if not j in DB.keys():
                missedkey = j
            if not j in DBCounter.keys():
                DBCounter[j] = 1
            else:
                DBCounter[j] += 1
    if missedkey != None:
        DB[missedkey] = []
    nodes = {}              # Store the imbalanced nodes.
    for key, val in DBCounter.items():
        if val > len(DB[key]):      # More incoming than outgoing - end here
            nodes['end'] = key
        elif val < len(DB[key]):    # More outgoing than incoming - start here.
            nodes['start'] = key
    DB[nodes['end']].append(nodes['start']) # Connect end to start to form Eulerian Cycle.
    cycle = EulerianCycle(DB)               # Circle up the Eulerian cycle.
    start = cycle.index(nodes['start'])      # Find the chopped point where one arrow is arbiturarily added.
    path = cycle[start:] + cycle[:start]    # Get the final path.
    return path



# Print out the Euler cycle or path.
def PrintEuler(DB, cycle = False):
    if cycle:
        path = EulerianCycle(DB)
        for i in path:
            print(i, '->', sep = '', end = '')
        print(path[0])
    else:
        path = EulerianPath(DB)
        print(path[0], end = '')
        for i in path[1:]:
            print('->', i, sep = '', end = '')
        


# Genome assembly using above functions.
def StrReconst(patterns):
    DB = GS1.DeBruijnP(patterns)
    path = EulerianPath(DB)
    genome = path[0]
    for i in path[1:]:
        genome += i[-1]
    return genome



# Solve k-Universal string problem.
def kUniStr(k):
    dim = 2**k
    coll = []
    model = '{0:0' + str(k) + 'b}'
    for i in range(dim):
        coll.append(model.format(i))
    DB = GS1.DeBruijnP(coll)
    cycle = EulerianCycle(DB)
    pattern = ''
    for i in cycle:
        pattern += i[0]
    return pattern



# Get (k,d)-mer of a genome sequence.
# e.g., the 1st (3,2)-mer of genome "GATACTAGACCG" is "GAT|TAG"
def kdmer(k, d, genome):
    l = k*2 + d
    n = len(genome)
    if n < l:
        print('Invalid (k,d)-mer!')
        return None
    coll = []
    for i in range(n-l+1):
        k1 = genome[i:(i+k)]
        k2 = genome[(i+k+d):(i+l)]
        coll.append(k1 + '|' + k2)
    coll.sort()
    return coll



# DeBruijn graph from (k,d)-mers pattern set, rather than single k-mer.
def PairedDeBruijn(k, patterns):
    DB = dict()
    if type(patterns[0]) is str:
        for i in range(len(patterns)):
            patterns[i] = patterns[i].split('|')    # Reformat 'k1|k2' into ['k1', 'k2'].
    for pattern in patterns:
        pre = pattern[0][:-1] + '|' + pattern[1][:-1]
        suf = pattern[0][1:] + '|' + pattern[1][1:]
        if not pre in DB.keys():
            DB[pre] = [suf]
        else:
            DB[pre].append(suf)
    return DB



# Genome assembly using above functions.
def StrReconstFromPair(k, d, patterns):
    DB = PairedDeBruijn(k, patterns)
    cache = cp.deepcopy(DB)
    while True:
        DB = cp.deepcopy(cache)
        path = EulerianPath(DB)
        pre = path[0][:(k-1)]
        suf = path[0][(1-k):]
        for pair in path[1:]:
            pre += pair[k-2]
            suf += pair[-1]
        if pre[(k+d):] == suf[:-(k+d)]:
            break
    genome = pre + suf[-(k+d):]
    return genome



# Find non-branched path collections from De Bruijn graph.
def MaxNonBranchingPath(graph):
    DB = ReadDB(graph)
    DBCounter = dict()     # This dict stores numbers of in-paths of items.
    for key, val in DB.items():
        if not key in DBCounter.keys():
            DBCounter[key] = [0, len(val)]
        for j in val:
            if not j in DB.keys():
                DBCounter[j] = [1, 0]
            if not j in DBCounter.keys():
                DBCounter[j] = [1, len(DB[j])]
            else:
                DBCounter[j][0] += 1
    paths = []
    cycle = []
    for key, val in DBCounter.items():
        if val[0] != 1 or val[1] != 1:  # Not 1-in-1-out node.
            if val[1] > 0:              # Can be starting node if outdegree is no zero.
                for nxtnode in DB[key]: # Iterate thru all outgoing paths until find non 1-in-1-out node.
                    path = [key]
                    cursor = nxtnode
                    while DBCounter[cursor][0] == 1 and DBCounter[cursor][1] == 1:
                        if cursor in path:
                            break
                        path.append(cursor)
                        cursor = DB[cursor][0]
                    path.append(cursor)
                    paths.append(path)
        else:
            path = [key]
            nxtnode = DB[key][0]
            while DBCounter[nxtnode][0] == 1 and DBCounter[nxtnode][1] == 1:
                path.append(nxtnode)
                nxtnode = DB[nxtnode][0]
                if nxtnode == key:
                    path.append(nxtnode)
                    exist = False
                    for item in cycle:
                        if set(path) == set(item):
                            exist = True
                            break
                    if not exist:
                        cycle.append(path)
                    break
    paths.extend(cycle)
    return paths
    


# Plot the non-branched path results.
def PlotNonBranPath(graph):
    paths = MaxNonBranchingPath(graph)
    for path in paths:
        for i in path[:-1]:
            print(i, '->', end = ' ')
        print(path[-1])


   
# Construct contig collection from patterns.
def ContigPatterns(patterns):
    DB = GS1.DeBruijnP(patterns)
    MaxPaths = MaxNonBranchingPath(DB)
    contigs = []
    for path in MaxPaths:
        contig = path[0]
        for i in path[1:]:
            contig += i[-1]
        contigs.append(contig)
    return contigs
    
