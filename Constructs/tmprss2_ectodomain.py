from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("NC_000021.9")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NC_000021.9"], "NC_000021.9", 'MALNSGSPPA',  min_peptide_length=492)

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_000021.9"], 106, 492, "TMPRSS2_ECD")

notISeq = tools.create_seq_object_from_string("GCGGCCGC")
tools.deconstruct_dna_sequence(notISeq, "NotI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

linker_0_seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker_0_seq, "Linker_0", True)

seqPeptideSeq = tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SecretionSignal_mouseIgKappa", True)

linker_0_seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker_0_seq, "Linker_0", True)

linker_1_mclellan = tools.create_seq_object_from_string(
    '''
    GGATCCGGA
    '''
)

tools.deconstruct_dna_sequence(linker_1_mclellan, 'Linker_1', True)

apTagSeq = tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstruct_dna_sequence(apTagSeq, "APTag", True)

hrv3c_protease_clevage_site = tools.create_seq_object_from_string(
    '''
    CTGGAGGTGCTGTTCCAGGGCCCA
    '''
)

tools.deconstruct_dna_sequence(hrv3c_protease_clevage_site, 'hrv3c_protease_cleavage', True)

eight_x_his_tag = tools.create_seq_object_from_string(
    '''
    CATCACCACCATCACCACCATCAT
    '''
)

tools.deconstruct_dna_sequence(eight_x_his_tag, 'his_8', True)

linker_4_mclellan = tools.create_seq_object_from_string(
    '''
    AGCGCC
    '''
)

tools.deconstruct_dna_sequence(linker_4_mclellan, 'Linker_4', True)

linker_2_mclellan = tools.create_seq_object_from_string(
    '''
    AGATCC
    '''
)

tools.deconstruct_dna_sequence(linker_2_mclellan, 'Linker_2', True)

stop_3 = tools.create_seq_object_from_string(
    '''
    TGATAATGA
    '''
)

tools.deconstruct_dna_sequence(stop_3, 'stops', True)

tools.create_construct_from_deconstructed_sequences(['SecretionSignal_mouseIgKappa', 'Linker_0', 'APTag', "Linker_1", 'TMPRSS2_ECD', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops'], 'BAP_TMPRSS2ECD_HRV3CProtease_His8x')
tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'APTag', "Linker_1", 'TMPRSS2_ECD', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops', 'XbaI'], 'BAP_TMPRSS2ECD_HRV3CProtease_His8x')