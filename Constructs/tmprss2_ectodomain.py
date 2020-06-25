from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("NM_001135099.1")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NM_001135099.1"], "NM_001135099.1", 'MALNSGSPPA',  min_peptide_length=492)

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NM_001135099.1"], 106, 492, "TMPRSS2_ECD")

gstSeq = tools.create_seq_object_from_string("ATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTTTTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAATGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGATGGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAACATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTGGATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTTGATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAAACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGATGTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAAAAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCATGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAA")
tools.deconstruct_dna_sequence(gstSeq, "GST", True)

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

tools.deconstruct_dna_sequence(linker_3_mclellan, 'Linker_3', True)

tools.deconstruct_dna_sequence(stop_3, 'stops', True)

tmprss2_ecd = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'FLAG3x', 'Linker_0', 'APTag', "Linker_1", 'TMPRSS2_ECD', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']

tools.create_construct_from_deconstructed_sequences(['SecretionSignal_mouseIgKappa', 'FLAG3x', 'Linker_0', 'APTag', "Linker_1", 'TMPRSS2_ECD', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_3', 'his_8', 'stops'], 'Flag_3x_BAP_TMPRSS2ECD_HRV3CProtease_His8x')
tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal_mouseIgKappa', 'FLAG3x', 'Linker_0', 'APTag', "Linker_1", 'TMPRSS2_ECD', 'Linker_2', 'hrv3c_protease_cleavage', "Linker_3", "GST", 'Linker_4', 'his_8', 'stops', 'XbaI'], 'NotI_Flag3x_BAP_TMPRSS2ECD_HRV3CProtease_GST_His8x_XbaI')

check_sequence = tools.create_seq_object_from_string(
    """
    MALNSGSPPAIGPYYENHGYQPENPYPAQPTVVPTVYEVHPAQYYPSPVPQYAPRVLTQA
SNPVVCTQPKSPSGTVCTSKTKKALCITLTLGTFLVGAALAAGLLWKFMGSKCSNSGIEC
DSSGTCINPSNWCDGVSHCPGGEDENRCVRLYGPNFILQVYSSQRKSWHPVCQDDWNENY
GRAACRDMGYKNNFYSSQGIVDDSGSTSFMKLNTSAGNVDIYKKLYHSDACSSKAVVSLR
CIACGVNLNSSRQSRIVGGESALPGAWPWQVSLHVQNVHVCGGSIITPEWIVTAAHCVEK
PLNNPWHWTAFAGILRQSFMFYGAGYQVEKVISHPNYDSKTKNNDIALMKLQKPLTFNDL
VKPVCLPNPGMMLQPEQLCWISGWGATEEKGKTSEVLNAAKVLLIETQRCNSRYVYDNLI
TPAMICAGFLQGNVDSCQGDSGGPLVTSKNNIWWLIGDTSWGSGCAKAYRPGVYGNVMVF
TDWIYRQMRADG
    """
)

check_sequence_2_ecd = tools.create_seq_object_from_string(
    """
    WKFMGSKCSNSGIECDSSGTCINPSNWCDGVSHCPGGEDENRCVRLYGPNFILQVYSSQR
KSWHPVCQDDWNENYGRAACRDMGYKNNFYSSQGIVDDSGSTSFMKLNTSAGNVDIYKKL
YHSDACSSKAVVSLRCIACGVNLNSSRQSRIVGGESALPGAWPWQVSLHVQNVHVCGGSI
ITPEWIVTAAHCVEKPLNNPWHWTAFAGILRQSFMFYGAGYQVEKVISHPNYDSKTKNND
IALMKLQKPLTFNDLVKPVCLPNPGMMLQPEQLCWISGWGATEEKGKTSEVLNAAKVLLI
ETQRCNSRYVYDNLITPAMICAGFLQGNVDSCQGDSGGPLVTSKNNIWWLIGDTSWGSGC
AKAYRPGVYGNVMVFTDWIYRQMRADG
    """
)



test_result = tools.compare_peptide_construct_to_sequence(tools.all_constructs['NM_001135099.1'], check_sequence)
test_result_2 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['TMPRSS2_ECD'], check_sequence_2_ecd)

pass