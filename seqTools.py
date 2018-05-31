def concatenate(alignments, padding_length=0):

    '''
    Source: https://github.com/karolisr/krpy/blob/master/kralign.py
    Concatenate alignments based on the Seq ids; row order does not
    matter. If one alignment contains a Seq id that another one does
    not, gaps will be introduced in place of the missing Seq.

    Args:
        alignments: (tuple, list) Alignments to be concatenated.

        padding_length: Introduce this many gaps between concatenated
            alignments.
    '''

    from Bio import Alphabet
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    if not isinstance(alignments, (list, tuple)):
        raise ValueError('Argument must be a list or a tuple.')
    elif len(alignments) == 1:
        return alignments[0]
    if isinstance(alignments, tuple):
        alignments = list(alignments)
    aln1 = None
    aln2 = None
    if len(alignments) > 2:
        aln2 = alignments.pop()
        aln1 = concatenate(alignments=alignments,
                           padding_length=padding_length)
    elif len(alignments) == 2:
        aln1 = alignments[0]
        aln2 = alignments[1]
    if (not isinstance(aln1, MultipleSeqAlignment) or
            not isinstance(aln2, MultipleSeqAlignment)):
        raise ValueError(
            'Argument must inherit from Bio.Align.MultipleSeqAlignment.')
    alphabet = Alphabet._consensus_alphabet([aln1._alphabet, aln2._alphabet])
    aln1_dict = dict()
    aln2_dict = dict()
    for aln1_s in aln1:
        aln1_dict[aln1_s.id] = aln1_s
    for aln2_s in aln2:
        aln2_dict[aln2_s.id] = aln2_s
    aln1_length = aln1.get_alignment_length()
    aln2_length = aln2.get_alignment_length()
    aln1_gaps = SeqRecord(Seq('-' * aln1_length, alphabet))
    aln2_gaps = SeqRecord(Seq('-' * aln2_length, alphabet))
    padding = SeqRecord(Seq('-' * padding_length, alphabet))
    result_seq_list = list()
    for aln1_key in aln1_dict.keys():
        merged_Seq = None
        if aln1_key in aln2_dict:
            merged_Seq = aln1_dict[aln1_key] + padding + aln2_dict[aln1_key]
            merged_Seq.id = aln1_dict[aln1_key].id
            merged_Seq.name = ''
            merged_Seq.description = ''
            aln2_dict.pop(aln1_key)
        else:
            aln1_seq_record = aln1_dict[aln1_key]
            merged_Seq = aln1_seq_record + padding + aln2_gaps
            merged_Seq.id = aln1_seq_record.id
            merged_Seq.name = ''
            merged_Seq.description = ''
        result_seq_list.append(merged_Seq)
    for aln2_seq_record in aln2_dict.values():
        merged_Seq = aln1_gaps + padding + aln2_seq_record
        merged_Seq.id = aln2_seq_record.id
        merged_Seq.name = ''
        merged_Seq.description = ''
        result_seq_list.append(merged_Seq)
    result_alignment = MultipleSeqAlignment(result_seq_list, alphabet)
    result_alignment.sort()
    return result_alignment




def genPartition(alignments, filename):
    # Generates a partition file based on the sequence lengths of multiple sequence alignments in a list
    # Args:
    #   alignments = List of AlignIO objects
    #   filename = Directory and filename of partition file
    # Returns
    #   Text file, partitions will be in the same sequence of alignments and is in the format for raxml
    
    import numpy as np
    
    genelen=[len(alignment[0]) for alignment in alignments]
    ends=np.cumsum(genelen)
    starts=[1]
    starts2=ends[0:len(ends)-1]+1 # keep in array form to do operations
    starts.extend(starts2.tolist())
    ends=ends.tolist() # convert back to list

    with open(filename, 'wb') as file:
        for i in range(len(starts)):
            num = i+1
            # if i == len(starts):
            #   file.write('DNA, gene' + str(num) + " = " + str(starts[i]) + '-' + str(ends[i]))
            # else:
            file.write('DNA, gene' + str(num) + " = " + str(starts[i]) + '-' + str(ends[i]) + '\n')