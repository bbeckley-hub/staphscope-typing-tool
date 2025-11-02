#!/usr/bin/env python3

import spaTyper.utils

####################################################
def check_primers(qDict):
    """
    Progress through a set of primers looking for an enriched sequence.
    
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    """

    seq_list = []
    for i in qDict.keys():
        enriched_seqs = enrichSeq(qDict[i].upper(), 'TAAAGACGATCCTTCGGTGAG', 'CAGCAGTAGTGCCGTTTGCTT')
        seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict.keys():
            enriched_seqs = enrichSeq(qDict[i].upper(), 'AGACGATCCTTCGGTGAGC', 'GCTTTTGCAATGTCATTTACTG')
            seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict.keys():
            enriched_seqs = enrichSeq(qDict[i].upper(), 'ATAGCGTGATTTTGCGGTT', 'CTAAATATAAATAATGTTGTCACTTGGA')
            seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict.keys():
            enriched_seqs = enrichSeq(qDict[i].upper(), 'CAACGCAATGGTTTCATCCA', 'GCTTTTGCAATGTCATTTACTG')
            seq_list += enriched_seqs
    if seq_list == []:
        return ['no enriched sequence.']
    if len(seq_list) > 1:
        sys.stderr.write(' more than one enriched sequence in ' + infile + '\n')
        
    return(seq_list)

####################################################
def enrichSeq(seq, fortemp, revtemp):
    """
    Simulates pcr on a sequence.
    
    This function takes sequences, forward pcr template and reverse pcr template.
    
    :param seq: Nucleotide sequence.
    :param fortemp: Forward primer sequence
    :param revtemp: Reverse primer sequence
    
    :type seq: string
    :type fortemp: string
    :type revtemp: string
    
    :returns: List of enriched sequences.
        
    .. attention:: Be aware of Copyright

        The code implemented here was retrieved and modified from spa_typing (https://github.com/mjsull/spa_typing)

        Give him credit accordingly.
    
    """
    index = 0
    found = None
    out1 = []
    out_list = []
    while found != -1:
        found = seq.find(fortemp, index)
        index = found + 1
        if found != -1:
            out1.append(found)
    index = 0
    found = None
    out2 = []
    revtempr = spaTyper.utils.revseq(revtemp)
    while found != -1:
        found = seq.find(revtempr, index)
        index = found + 1
        if found != -1:
            out2.append(found)
    for j in out1:
        for k in out2:
            if k - j > 50 and k - j < 5000:
                enriched_seq = seq[j:k+len(revtemp)]
                out_list.append(enriched_seq)
    index = 0
    found = None
    out1 = []
    fortempr = spaTyper.utils.revseq(fortemp)
    while found != -1:
        found = seq.find(fortempr, index)
        index = found + 1
        if found != -1:
            out1.append(found)
    index = 0
    found = None
    out2 = []
    while found != -1:
        found = seq.find(revtemp, index)
        index = found + 1
        if found != -1:
            out2.append(found)
    for j in out2:
        for k in out1:
            if k - j > 50 and k - j < 5000:
                enriched_seq = seq[j:k+len(fortemp)]
                out_list.append(spaTyper.utils.revseq(enriched_seq))
    return out_list
