# -*- coding: utf-8 -*-
"""
Genome Sequencing 3

Created on Thu Feb  2 14:03:40 2017

@author: GY
"""

import GS2



# Get compliment nucleotide.
def compnuc(char):
    if char=='A' or char=='a':
        return 'T'
    elif char=='T' or char=='t':
        return 'A'
    elif char=='C' or char=='c':
        return 'G'
    elif char=='G' or char=='g':
        return 'C'
    else:
        #print('Error: wrong neucliotide letter!')
        return ''


    
# Return compliment strand of a given DNA sequence.
def compstr(text):
    comp=''
    for i in text:
        comp=compnuc(i)+comp
    return comp




# Masses of single aminoacid
AminoAcidMass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113,
'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}

Inv_AminoAcidMass = {}
for k, v in AminoAcidMass.items():
    if not v in Inv_AminoAcidMass.keys():
        Inv_AminoAcidMass[v] = [k]
    else:
        Inv_AminoAcidMass[v].append(k)


# Codon dictionary
codondict = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T',
 'ACG': 'T', 'ACU': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I',
 'AUG': 'M', 'AUU': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P',
 'CCG': 'P', 'CCU': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L',
 'CUG': 'L', 'CUU': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D', 'GCA': 'A', 'GCC': 'A',
 'GCG': 'A', 'GCU': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V',
 'GUG': 'V', 'GUU': 'V', 'UAA': '', 'UAC': 'Y', 'UAG': '', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S',
 'UCG': 'S', 'UCU': 'S', 'UGA': '', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F',
 'UUG': 'L', 'UUU': 'F'}

Inv_codondict = {}
for k, v in codondict.items():
    if not v in Inv_codondict.keys():
        Inv_codondict[v] = [k]
    else:
        Inv_codondict[v].append(k)



# Translate RNA sequence to protein sequence.
def RNAtoPiptide(RNA):
    n = len(RNA)
    pep = ''
    idx = 0
    while idx < n//3:
        codon = RNA[(idx*3):(idx*3+3)]
        if not codondict[codon]:
            break
        pep += codondict[codon]
        idx += 1
    return pep



# Find potential peptide encoding sites from genome sequence.
def PeptideEncodingPats(seq, pep, reverse = True):
    if 'T' in seq:          # Change DNA to RNA.
        seq = seq.replace('T', 'U')
    seqset = ['']
    if reverse:
        compset = ['']
    for aa in pep:                          # Scan thru amino acids in given peptide.
        cache = []
        if reverse:
            compcache = []
        for codon in Inv_codondict[aa]:     # Scan thru each codon for the given acid.
            for codons in seqset:          # Scan thru each pattern of codons adding the new codon.
                newcodons = codons + codon
                if newcodons in seq:        # Check the newly formed codons. If in, the add into list.
                    cache.append(newcodons)
            if reverse:
                for compcodons in compset:
                    newcodons = compcodons + codon
                    if compstr(newcodons.replace('U', 'T')).replace('T', 'U') in seq:
                        compcache.append(newcodons)
        seqset = cache
        if reverse:
            compset = compcache
    seqset = set(seqset)
    if reverse:
        compset = set(compset)
    seqlist = []
    k = len(pep) * 3
    for i in range(len(seq) - k + 1):
        for code in seqset:
            if code == seq[i:(i+k)]:
                seqlist.append(code.replace('U', 'T'))
                break
        for code in compset:
            rawcode = compstr(code.replace('U', 'T'))
            if rawcode.replace('T', 'U') == seq[i:(i+k)]:
                seqlist.append(rawcode)
    return seqlist



# Find theoretical spectrum of linear peptide.
def LinearSpectrum(pep, cyclic = False):
    spec = [0]
    n = len(pep)
    for i in range(n):          # Iterate by shifting starting point.
        mass = 0
        for j in range(i, n):   # Collect masses of subpeptides starting at ith position.
            mass += AminoAcidMass[pep[j]]
            spec.append(mass)
    spec.sort()
    if cyclic:
        molmass = spec[-1]      # Find molecular mass (last one in sorted list).
        exclude = LinearSpectrum(pep[1:-1]) # Find all substrings excluding head and tail.
        exclude.pop(0)          # excluding mass 0 - it's the peptide mass counted above.
        for excmass in exclude:
            spec.append(molmass - excmass)  # These masses are sub-peptides including first and last amino acids.
        spec.sort()
    return spec



# Find possible peptides from spectra.
def CycloPeptideSequencing(spec):
    if type(spec[0]) is str:        # Switch strings to numbers.
        specint = []
        for i in spec:
            specint.append(int(i))
        spec = specint
    peps = ['']                     # Initialize peptide collections
    spec.sort()
    def inlist(subspec, spec):      # Define function to compare subspectrum and spectrum.
        for i in subspec:
            if not i in spec:
                return False
        return True
    while True:
        cache = peps.copy()
        peps = []
        for pep in cache:           # Scan consistent peptides after adding each amino acid.
            for aa in Inv_AminoAcidMass.values():
                newpep = pep + aa[0]
                newspec = LinearSpectrum(newpep)
                if inlist(newspec, spec):
                    peps.append(newpep)
        if max(LinearSpectrum(peps[0])) == spec[-1]:
            break
    spectra = []
    for pep in peps:                # Transfer peptides to spectra.
        newspec = []
        for aa in pep:
            newspec.append(AminoAcidMass[aa])
        spectra.append(newspec)
    return spectra
                
























