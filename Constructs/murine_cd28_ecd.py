from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")

sequence_identifier = "NM_007642.4"

tools.import_sequence_by_ncbi_identifier(sequence_identifier)
tools.deconstruct_imported_cdna_sequence(tools.all_sequences[sequence_identifier], sequence_identifier, 218)
tools.deconstruct_imported_orf_sequence(tools.all_sequences[sequence_identifier], sequence_identifier, 'MTLRLLFLAL',  min_peptide_length=218)
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs[sequence_identifier], 20, 150, "CD28_Extracellular")
dnaSeq = tools.return_dna_sequence_from_deconstructed_list(tools.all_deconstructed_sequences["CD28_Extracellular"]['deconstructedList'])

nhe1Seq = tools.create_seq_object_from_string("GCTAGC")
tools.deconstruct_dna_sequence(nhe1Seq, "NheI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

ecoR1Seq = tools.create_seq_object_from_string("GAATTC")
tools.deconstruct_dna_sequence(ecoR1Seq, "EcoRI", False)

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

flag3x = tools.create_seq_object_from_string(
    '''
    gactacaaagaccatgacggtgattataaagatcatgacatcgattacaaggatgacgatgacaag
    '''
)

tools.deconstruct_dna_sequence(flag3x, 'FLAG3x', True)

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

linker_3_mclellan = tools.create_seq_object_from_string(
    '''
    GGC
    '''
)

linker_4_mclellan = tools.create_seq_object_from_string(
    '''
    AGCGCC
    '''
)

tools.deconstruct_dna_sequence(linker_3_mclellan, 'Linker_3', True)

# Specific linker before TST tag
tools.deconstruct_dna_sequence(linker_4_mclellan, 'Linker_4', True)

tools.deconstruct_dna_sequence(stop_3, 'stops', True)

kozak_seq_non_coding = tools.create_seq_object_from_string("gccAcc")
tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakNonCoding", False)

linker_5_seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker_5_seq, "Linker_5", True)

twin_strept_2_tag = tools.create_seq_object_from_string(
    '''
    TGGTCCCACCCCCAGTTCGAGAAGGGCGGCGGTAGTGGAGGGGGCGGATCTGGCGGCTCAGCTTGGAGCCACCCCCAGTTCGAAAAG
    '''
)

tools.deconstruct_dna_sequence(twin_strept_2_tag, "TST", True)

cd28_construct = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'CD28_Extracellular', "Linker_1", 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']

cd28_construct_w_re = cd28_construct.copy()
cd28_construct_w_re.insert(0, "NheI")
cd28_construct_w_re.append("XbaI")

tools.create_construct_from_deconstructed_sequences(cd28_construct, 'secretion_tst_his8_hrv3c_AP_CD28_extracellular')
tools.create_construct_from_deconstructed_sequences(cd28_construct_w_re, 'nheI_secretion_CD28_extracellular_AP_hrv3c_his8_tst_xbaI')

# cd40l_construct = ['Linker_5', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']

pass