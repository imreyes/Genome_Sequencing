# -*- coding: utf-8 -*-
"""
Genome Sequencing 4

Further optimization for practical peptide sequencing and reconstruction
from mass patterns

Created on Thu Feb  2 14:03:40 2017

@author: GY
"""

import GS3
import pandas as pd
import copy



# Score the given cyclo peptide sequence based on given mass pattern.
# Match a pair of same mass number, and give 1 pt. Return summed scores.
def pepScore(pep, spec, cyclic = False):
    RealSpec = spec.copy()
    RealSpec.sort()
    TheoSpec = GS3.LinearSpectrum(pep, cyclic = cyclic)
    TheoSpec.sort()
    score = 0
    while RealSpec and TheoSpec:
        R = RealSpec[0]
        T = TheoSpec[0]
        if R == T:
            RealSpec.pop(0)
            TheoSpec.pop(0)
            score += 1
        elif R < T:
            RealSpec.pop(0)
        else:
            TheoSpec.pop(0)
    return score







# Find best-matched peptide sequence with highest score.
# Applying Leader Board mechanism.
def LeaderboardCyclopeptideSequence(spec, N, cyclic = False):
    spec = list(map(int, spec.split()))
    # First establish peptide mass table with non-redundant masses.
    NonRedMass = {}         # Build peptide mass library - ignore multiplicity of amino acids with same MW.
    for key, value in GS3.Inv_AminoAcidMass.items():
        NonRedMass[value[0]] = key
    Leaderboard = pd.Series()
    for key in NonRedMass.keys():
        Leaderboard[key] = pepScore(key, spec, cyclic = cyclic)
    limit = max(spec)       # Parent mass - all candidates must have this total mass.
    while not Leaderboard.empty:
        Newboard = pd.Series()      # Build new board.
        for pep, score in Leaderboard.items():
            for aa in NonRedMass:
                newpep = pep + aa      # add one amino acid to the right.
                newmass = GS3.pepMass(newpep)
                if newmass < limit - 56 or newmass == limit:
                    newscore = pepScore(newpep, spec, cyclic = cyclic)
                    Newboard[newpep] = newscore     # Create new row for Leaderboard.
        if Newboard.empty:      # Empty Newboard means no update, and the while loop is break.
            break
        Leaderboard = Leaderboard.append(Newboard)      # Add Newboard's entries to Leaderboard.
        LBvol = len(Leaderboard)
        if len(Leaderboard) <= N:  # Volume of Leaderboard is kept if smaller than cut volume.
            continue
        Leaderboard.sort_values(ascending = False, inplace = True)  # Sort by score.
        cutscore = Leaderboard[N-1]     # Cut is set to score of N th peptide in board.
        start = N
        while Leaderboard[-1] != cutscore:
            LBvol = len(Leaderboard)
            cut = (LBvol - 1 + start) // 2
            if Leaderboard[cut] < cutscore:
                Leaderboard = Leaderboard[:cut]
            elif LBvol - cut <= 2:
                Leaderboard = Leaderboard[:-1]
            else:
                start = cut
        print(Leaderboard[0:1], end = ' ')
    for pep in Leaderboard.index:
        if GS3.pepMass(pep) != limit:
            Leaderboard.drop(pep, inplace = True)
    LeaderPep = Leaderboard.index[0]
    #print(LeaderPep)
    #print(type(LeaderPep))
    print(str(GS3.AminoAcidMass[LeaderPep[0]]), end = '')
    for aa in LeaderPep[1:]:
        print('-' + str(GS3.AminoAcidMass[aa]), end = '')
    return LeaderPep



# Calculate spectral convolution - all mass differences.
def SpectralConvolution(spec, M = 0):
    l = len(spec)
    spec.sort()
    conv = pd.Series()      # Use series to store leaderboard of convolution.
    for i in range(l-1):
        for j in range(i+1, l):
            diff = spec[j] - spec[i]
            if diff >= 57 and diff <= 200:      # Mass difference must be in range of amino acid masses.
                if diff in conv.keys():
                    conv.loc[diff] += 1
                else:
                    conv.loc[diff] = 1
    conv.sort_values(ascending = False, inplace = True)
    if M:
        cutscore = conv.iloc[M-1]     # Cut is set to score of N th mass in board.
        start = M
        while conv.iloc[-1] != cutscore:
            LBvol = len(conv)
            cut = (LBvol - 1 + start) // 2
            if conv.iloc[cut] < cutscore:
                conv = conv[:cut]
            elif LBvol - cut <= 2:
                conv = conv[:-1]
            else:
                start = cut
    return conv



# Find best-matched peptide sequence with highest score.
# Applying Leader Board mechanism.
def ConvolutionCyclopeptideSequence(spec, M, N, cyclic = False):
    if type(spec) is str:
        spec = list(map(int, spec.split()))
    spec.sort()
    # First establish peptide mass table with non-redundant masses.
    MassLib = SpectralConvolution(spec, M).index
    Leaderboard = []
    limit = max(spec)       # Parent mass - all candidates must have this total mass.
    for mass in MassLib:    # Initialize Leaderboard by adding 1 amino acid mass.
        resid = spec.copy()
        if 0 in resid:
            resid.remove(0)
        if mass in resid:    # if the mass of amino acid is in spec, then score is 2 (0, mass); otherwise 1 (0).
            resid.remove(mass)
            Leaderboard.append([1, [mass], resid])
        else:
            Leaderboard.append([0, [mass], resid])
    scorerec = []
    while Leaderboard:
        Newboard = []      # Build new board.
        for item in Leaderboard:
            for mass in MassLib:
                newitem = copy.deepcopy(item)
                newitem[1].append(mass)      # add one amino acid to the right.
                newmass = sum(newitem[1])
                if newmass < limit - 56 or newmass == limit:
                    for i in range(1, len(newitem[1])): # Edit newitem:
                        newpat = sum(newitem[1][-i:])   # New mass pattern
                        if newpat in newitem[2]:        # When found in residual spectrum
                            newitem[0] += 1             # Score add 1
                            newitem[2].remove(newpat)   # Mass removed from residual spectrum
                    Newboard.append(newitem)
        if not Newboard:      # Empty Newboard means no update, and the while loop is break.
            break
        Leaderboard.extend(Newboard)      # Add Newboard's entries to Leaderboard.
        LBvol = len(Leaderboard)
        if len(Leaderboard) <= N:  # Volume of Leaderboard is kept if smaller than cut volume.
            continue
        Leaderboard.sort(reverse = True)         # Sort by score.
        scorerec.append(Leaderboard[0][0])
        cutscore = Leaderboard[N-1][0]     # Cut is set to score of N th peptide in board.
        start = N
        while Leaderboard[-1][0] != cutscore:
            LBvol = len(Leaderboard)
            cut = (LBvol - 1 + start) // 2
            if Leaderboard[cut][0] < cutscore:
                Leaderboard = Leaderboard[:cut]
            elif LBvol - cut <= 2:
                Leaderboard = Leaderboard[:-1]
            else:
                start = cut
        if len(scorerec) > limit/200:
            if len(scorerec) >= 5:
                if scorerec[-1] == scorerec[-5]:
                    break
    for item in Leaderboard:
        if sum(item[1]) != limit:
            Leaderboard.remove(item)
    LeaderMasses = Leaderboard[0][1]
    print(LeaderMasses[0], end = '')
    for mass in LeaderMasses[1:]:
        print('-', mass, sep = '', end = '')
    return LeaderMasses







