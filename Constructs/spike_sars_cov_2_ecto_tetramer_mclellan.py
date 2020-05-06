from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")

# import spike sequence
tools.import_sequence_by_ncbi_identifier("NC_045512.2")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NC_045512.2"], "NC_045512.2", 'MFVFLVLLPL',  min_peptide_length=1273)

tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP', 986, 'K', 'P')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP_987_VP', 987, 'V', 'P')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP_987_VP_682_RG', 682, 'R', 'G')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP_987_VP_682_RG_683_RS', 683, 'R', 'S')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS', 685, 'R', 'S')

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 13, 1273, "Spike_ecto_domain_mclellan")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 13, 685, "Spike_S1_ecto_domain_mclellan")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 816, 1273, "Spike_S2_ecto_domain_mclellan")

notISeq = tools.create_seq_object_from_string("GCGGCCGC")
tools.deconstruct_dna_sequence(notISeq, "NotI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

linker_0_seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker_0_seq, "Linker_0", True)

seqPeptideSeq = tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SecretionSignal_mouseIgKappa", True)

sars_cov2_ecto_mclellan = tools.create_seq_object_from_string(
    '''
    ATGTTCGTGTTCCTGGTGCTCCTGCCTCTGGTGAGCAGCCAGTGCGTGAACCTGACCACCCGAACCCAGCTCC
    CACCAGCCTACACCAACAGCTTTACACGGGGCGTGTACTACCCTGACAAGGTGTTCAGATCTAGCGTCCTGCA
    CAGCACTCAGGACCTCTTCCTGCCGTTCTTCAGCAACGTGACATGGTTCCACGCCATCCACGTGAGCGGCACA
    AACGGAACCAAGCGGTTTGATAACCCCGTCCTGCCATTCAATGATGGAGTTTACTTCGCCAGTACCGAGAAGA
    GTAACATCATCCGGGGCTGGATCTTCGGCACCACCCTGGATAGCAAAACACAGAGCCTCCTGATCGTGAACAA
    TGCCACGAACGTCGTGATCAAGGTGTGCGAGTTCCAGTTTTGCAATGATCCTTTCCTGGGTGTGTACTACCAC
    AAGAACAACAAGAGCTGGATGGAAAGCGAGTTCAGAGTCTACAGCAGCGCCAACAACTGCACATTCGAGTACG
    TCTCTCAGCCTTTTCTGATGGACCTTGAGGGGAAACAAGGCAACTTCAAGAACCTGAGAGAATTCGTGTTCAA
    GAACATCGACGGCTACTTCAAAATCTACTCCAAGCACACACCCATCAACCTGGTCCGGGACCTCCCTCAGGGC
    TTCAGCGCCCTGGAACCCCTGGTCGACCTGCCCATAGGCATCAACATAACGCGGTTCCAAACCCTGCTGGCCC
    TGCATAGATCCTACCTGACTCCTGGCGACAGCAGCAGCGGATGGACCGCCGGAGCTGCAGCCTACTATGTGGG
    CTACCTGCAACCTAGAACCTTCCTGCTGAAGTACAACGAGAACGGCACAATCACAGACGCCGTCGACTGCGCC
    CTGGACCCTCTCTCTGAGACAAAGTGCACCCTGAAGTCCTTCACCGTGGAAAAGGGCATCTACCAGACCAGCA
    ACTTCCGGGTGCAGCCTACAGAGAGCATCGTGCGATTTCCAAACATTACCAACCTCTGCCCCTTCGGCGAGGT
    GTTTAACGCCACAAGATTTGCCTCCGTTTACGCCTGGAATAGAAAGAGAATCAGCAATTGTGTGGCCGACTAC
    TCCGTGCTGTATAACAGCGCCTCTTTCAGCACCTTCAAGTGCTACGGCGTTTCCCCAACAAAGCTGAATGACC
    TGTGCTTCACCAACGTGTACGCCGACTCCTTCGTAATTAGAGGCGATGAGGTGCGGCAGATCGCACCAGGCCA
    GACCGGTAAGATCGCTGACTACAACTATAAGCTGCCTGATGATTTTACAGGCTGCGTGATCGCCTGGAACTCT
    AACAACCTGGATAGCAAGGTGGGCGGCAACTACAACTACCTGTACCGGCTGTTTCGCAAGTCTAACCTGAAAC
    CTTTCGAGAGAGACATCTCCACAGAGATCTACCAGGCCGGTTCTACACCTTGTAACGGGGTGGAAGGCTTCAA
    CTGTTACTTCCCTCTGCAAAGCTACGGCTTCCAGCCTACCAATGGAGTCGGCTACCAGCCATACCGGGTGGTC
    GTGCTGTCCTTCGAGTTACTCCACGCCCCCGCCACCGTCTGCGGTCCTAAGAAGTCCACCAATCTGGTTAAGA
    ACAAATGCGTGAACTTCAACTTCAACGGCCTGACCGGGACCGGCGTGCTGACCGAAAGCAACAAAAAGTTCCT
    CCCCTTCCAGCAGTTCGGCCGTGATATCGCTGACACCACAGATGCCGTCAGAGATCCACAGACCCTGGAAATC
    CTGGATATTACACCCTGCTCCTTCGGAGGAGTTTCTGTGATCACCCCCGGGACCAATACCAGCAACCAGGTGG
    CTGTGCTGTACCAAGATGTTAACTGCACCGAGGTTCCTGTGGCCATCCACGCCGATCAGCTGACACCTACTTG
    GAGAGTGTACTCCACTGGCTCCAATGTGTTCCAGACCAGGGCCGGATGTCTGATCGGCGCCGAGCACGTGAAT
    AACAGTTACGAGTGCGACATCCCTATCGGCGCCGGCATCTGTGCCAGCTACCAGACCCAGACAAACAGCCCTG
    GGTCTGCTTCCTCTGTAGCTAGCCAGAGCATCATCGCCTACACCATGAGCCTGGGCGCAGAGAACAGCGTGGC
    CTATTCCAACAACTCTATCGCCATTCCCACCAACTTTACAATTAGCGTCACAACAGAGATCCTGCCCGTGAGC
    ATGACCAAGACCAGCGTGGACTGTACAATGTACATCTGTGGCGACAGCACTGAATGCAGCAACCTGCTGCTGC
    AATACGGCTCCTTTTGCACCCAACTGAACCGGGCGCTGACCGGAATCGCCGTGGAACAGGACAAAAATACCCA
    GGAGGTGTTCGCCCAAGTGAAGCAGATCTACAAGACCCCACCTATCAAGGACTTCGGCGGCTTTAACTTTAGC
    CAGATTCTCCCTGATCCTTCTAAGCCTAGCAAGCGGAGCTTTATCGAGGATCTGCTGTTCAACAAGGTCACCC
    TGGCCGATGCCGGCTTTATCAAACAGTATGGCGATTGCCTGGGCGACATAGCCGCCAGAGATCTGATCTGCGC
    CCAGAAATTCAACGGCCTGACAGTTCTCCCACCTCTGCTGACCGACGAGATGATCGCTCAGTACACCTCTGCC
    CTGCTGGCTGGCACCATCACATCTGGGTGGACATTTGGCGCCGGCGCCGCCCTGCAGATCCCCTTTGCCATGC
    AGATGGCCTATAGATTCAACGGAATCGGCGTGACCCAGAACGTGCTGTATGAAAACCAGAAGCTGATCGCTAA
    CCAGTTCAATTCTGCCATCGGCAAGATCCAGGACTCCCTCTCCTCTACCGCCAGCGCCCTGGGCAAACTGCAG
    GACGTGGTGAATCAGAACGCCCAAGCCCTGAACACCCTGGTGAAGCAGCTCAGCAGCAATTTTGGCGCCATCA
    GCTCTGTGCTGAACGATATCCTGTCTAGACTGGACCCTCCAGAAGCCGAAGTCCAGATCGATAGACTGATCAC
    AGGCAGACTGCAGTCCCTGCAAACCTACGTGACCCAACAGCTGATCAGGGCCGCTGAAATAAGAGCCAGCGCC
    AATCTCGCCGCTACCAAGATGTCCGAGTGTGTGCTGGGACAGTCTAAACGCGTTGACTTCTGCGGCAAAGGCT
    ATCACCTGATGAGCTTCCCCCAGAGCGCGCCGCACGGCGTGGTGTTCCTGCATGTGACATACGTGCCTGCCCA
    AGAGAAGAATTTCACAACCGCCCCTGCCATCTGCCACGACGGCAAGGCCCACTTCCCAAGAGAGGGCGTTTTC
    GTTTCCAATGGCACACACTGGTTCGTGACACAAAGAAACTTCTACGAACCCCAGATTATCACCACCGACAACA
    CCTTCGTGAGTGGCAATTGTGACGTGGTCATCGGAATCGTGAACAACACAGTGTACGACCCTCTGCAACCTGA
    GCTGGACTCTTTTAAGGAAGAGCTGGACAAGTACTTTAAAAACCACACCAGCCCCGATGTGGACCTGGGCGAC
    ATCAGTGGCATTAACGCCAGCGTGGTGAACATCCAAAAGGAAATCGACAGACTGAACGAGGTGGCCAAGAACC
    TGAACGAGTCCCTGATCGACCTGCAGGAGCTCGGCAAATACGAGCAG
    '''
)

linker_1_mclellan = tools.create_seq_object_from_string(
    '''
    GGATCCGGA
    '''
)

t4_fibritin_trimerization_motif = tools.create_seq_object_from_string(
    '''
    TACATCCCCGAGGCCCCCAGAGATGGCCAGGCCTACGTGCGGAAGGACGGCGAGTGGGTACTGCTGAGCACATTCCTGGGC
    '''
)

linker_2_mclellan = tools.create_seq_object_from_string(
    '''
    AGATCC
    '''
)

hrv3c_protease_clevage_site = tools.create_seq_object_from_string(
    '''
    CTGGAGGTGCTGTTCCAGGGCCCA
    '''
)

linker_3_mclellan = tools.create_seq_object_from_string(
    '''
    GGC
    '''
)

eight_x_his_tag = tools.create_seq_object_from_string(
    '''
    CATCACCACCATCACCACCATCAT
    '''
)

linker_4_mclellan = tools.create_seq_object_from_string(
    '''
    AGCGCC
    '''
)

twin_strept_2_tag = tools.create_seq_object_from_string(
    '''
    TGGTCCCACCCCCAGTTCGAGAAGGGCGGCGGTAGTGGAGGGGGCGGATCTGGCGGCTCAGCTTGGAGCCACCCCCAGTTCGAAAAG
    '''
)

stop_3 = tools.create_seq_object_from_string(
    '''
    TGATAATGA
    '''
)

apTagSeq = tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstruct_dna_sequence(apTagSeq, "APTag", True)

tools.deconstruct_dna_sequence(sars_cov2_ecto_mclellan, "SarsCov2EctoMclellan", True)
tools.deconstruct_dna_sequence(stop_3, 'end_to_hind_3', True)
tools.deconstruct_dna_sequence(linker_1_mclellan, 'Linker_1', True)
tools.deconstruct_dna_sequence(t4_fibritin_trimerization_motif, 'T4_fibritin_trimerization_motif', True)
tools.deconstruct_dna_sequence(linker_2_mclellan, 'Linker_2', True)
tools.deconstruct_dna_sequence(hrv3c_protease_clevage_site, 'hrv3c_protease_cleavage', True)
tools.deconstruct_dna_sequence(linker_4_mclellan, 'Linker_4', True)
tools.deconstruct_dna_sequence(eight_x_his_tag, 'his_8', True)
tools.deconstruct_dna_sequence(stop_3, 'stops', True)

tools.create_construct_from_deconstructed_sequences(['SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops'], 'spike_sars_cov_2_trimer_bap_hrv3c_his8')
tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops', 'XbaI'], 'NotI_spike_sars_cov_2_trimer_bap_hrv3c_his8_XbaI')

tools.create_construct_from_deconstructed_sequences(['SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_S1_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops'], 'spike_s1_sars_cov_2_trimer_bap_hrv3c_his8')
tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_S1_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops', 'XbaI'], 'NotI_spike_s1_sars_cov_2_trimer_bap_hrv3c_his8_XbaI')

tools.create_construct_from_deconstructed_sequences(['SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_S2_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops'], 'spike_s2_sars_cov_2_trimer_bap_hrv3c_his8')
tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_S2_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage',  'Linker_4', 'his_8', 'stops', 'XbaI'], 'NotI_spike_s2_sars_cov_2_trimer_bap_hrv3c_his8_XbaI')

pass