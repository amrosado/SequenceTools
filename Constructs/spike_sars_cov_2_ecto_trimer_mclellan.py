from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")

# import spike sequence
tools.import_sequence_by_ncbi_identifier("NC_045512.2")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NC_045512.2"], "NC_045512.2", 'MFVFLVLLPL',  min_peptide_length=1273)

tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_proline_986_KP', 986, 'K', 'P')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences["NC_045512.2"], 'Spike_D614G', 614, 'D', 'G')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences['Spike_proline_986_KP'], 'Spike_proline_986_KP_987_VP', 987, 'V', 'P')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences['Spike_proline_986_KP_987_VP'], 'Spike_proline_986_KP_987_VP_682_RG', 682, 'R', 'G')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences['Spike_proline_986_KP_987_VP_682_RG'], 'Spike_proline_986_KP_987_VP_682_RG_683_RS', 683, 'R', 'S')
tools.make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(tools.all_deconstructed_sequences['Spike_proline_986_KP_987_VP_682_RG_683_RS'], 'Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS', 685, 'R', 'S')

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_045512.2"], 13, 1273, "Spike_ecd")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_045512.2"], 13, 685, "Spike_ecd_S1")

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 13, 1273, "Spike_ecto_domain_mclellan")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 13, 685, "Spike_S1_ecto_domain_mclellan")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS"], 686, 1273, "Spike_S2_ecto_domain_mclellan")

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_045512.2"], 13, 816, "Spike_S1_to_S2prime")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_045512.2"], 816, 1273, "Spike_ecd_S2prime")
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NC_045512.2"], 686, 1273, "Spike_S2")


kozak_seq_non_coding = tools.create_seq_object_from_string("gccAcc")
tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakNonCoding", False)

kozak_seq_coding = tools.create_seq_object_from_string("ATGGGT")
tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakCoding", True)

notISeq = tools.create_seq_object_from_string("GCGGCCGC")
tools.deconstruct_dna_sequence(notISeq, "NotI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

linker_0_seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker_0_seq, "Linker_0", True)

linker_5_seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker_5_seq, "Linker_5", True)

# Secretion signal mouse Ig kappa chain
seqPeptideSeq = tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SS_MIgK", True)

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

flag_tag = tools.create_seq_object_from_string(
    '''
    GATTACAAGGATGACGACGATAAG 
    '''
)

tools.deconstruct_dna_sequence(flag_tag, 'FLAG', True)

myc3x = tools.create_seq_object_from_string(
    '''
    GAGCAGAAACTCATCTCTGAAGAAGATCTGGAACAAAAGTTGATTTCAGAAGAAGATCTGGAACAGAAGCTCATCTCTGAGGAAGATCTG
    '''
)

tools.deconstruct_dna_sequence(myc3x, 'MYC3x', True)

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
tools.deconstruct_dna_sequence(linker_3_mclellan, 'Linker_3', True)
tools.deconstruct_dna_sequence(eight_x_his_tag, 'his_8', True)
tools.deconstruct_dna_sequence(twin_strept_2_tag, "TST", True)
tools.deconstruct_dna_sequence(stop_3, 'stops', True)

spike_trimer = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecd', "Linker_1", 'T4_fibritin_trimerization_motif', 'Linker_5', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']
spike_trimer_mclellan = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecto_domain_mclellan', "Linker_1", 'T4_fibritin_trimerization_motif', 'Linker_5', 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']
spike_s1 = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecd_S1', "Linker_1", 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']
spike_s2_prime = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_ecd_S2prime', "Linker_1", 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']
spike_s1_to_S2prime = ['KozakNonCoding', 'SecretionSignal_mouseIgKappa', 'Linker_0', 'Spike_S1_to_S2prime', "Linker_1", 'APTag', 'Linker_2', 'hrv3c_protease_cleavage', 'Linker_3', 'his_8', 'Linker_4', 'TST', 'stops']

spike_trimer_with_re = spike_trimer.copy()
spike_trimer_mclellan_with_re = spike_trimer_mclellan.copy()
spike_s1_with_re = spike_s1.copy()
spike_s2_prime_with_re = spike_s2_prime.copy()

tools.create_construct_from_deconstructed_sequences(spike_trimer, 'SC2_Spike_ecd_trimer')
tools.create_construct_from_deconstructed_sequences(spike_trimer_mclellan, 'SC2_Spike_ecd_trimer_mclellan')
tools.create_construct_from_deconstructed_sequences(spike_s1, 'SC2_Spike_S1')
tools.create_construct_from_deconstructed_sequences(spike_s2_prime, 'SC2_Spike_S2prime')
tools.create_construct_from_deconstructed_sequences(spike_s1_to_S2prime, 'SC2_Spike_S1ToS2prime')

compare_1 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['SarsCov2EctoMclellan'], tools.all_deconstructed_sequences['Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS']['peptideSequence'])


check_s1_protein_seq = tools.create_seq_object_from_string(
    '''
    SQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSG
TNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCE
FQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREF
VFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPG
DSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGI
YQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSA
SFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGC
VIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL
QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVL
TESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLY
QDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
SYQTQTNSPRRAR
    '''
)

check_s2_peptide_seq = tools.create_seq_object_from_string(
    '''
    SVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGD
STECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQI
LPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLL
TDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIAN
QFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLD
KVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGK
GYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVT
QRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVD
LGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLI
AIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
    '''
)

check_s2prime_peptide_seq = tools.create_seq_object_from_string(
    '''
    SFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTS
ALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQ
DSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDR
LITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQS
APHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQII
TTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINAS
VVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLC
CMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
    '''
)

peptide_seq_s_ecto_dfu_2p = tools.create_seq_object_from_string(
    '''
    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQ
    '''
)

# test = tools.compare_peptide_construct_to_sequence(tools.all_constructs['Spike_proline_986_KP_987_VP_682_RG_683_RS_685_RS'], peptide_seq_s_ecto_dfu_2p)
# test_s1 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['Spike_S1_ecto_domain_mclellan'], check_s1_protein_seq)
# test_s2 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['Spike_S2'], check_s2_peptide_seq)
# test_s2prime = tools.compare_peptide_construct_to_sequence(tools.all_constructs['Spike_S2_prime'], check_s2prime_peptide_seq)

pass