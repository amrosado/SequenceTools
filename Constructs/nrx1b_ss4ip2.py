from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")

nhe1Seq = tools.create_seq_object_from_string("GCTAGC")
tools.deconstruct_dna_sequence(nhe1Seq, "NheI", False)

nrx1b = tools.create_seq_object_from_string("""
ATGTACCAGAGGATGCTTCGGTGCGGCGCCGATCTGGGATCGCCCGGGGGCGGCAGTGGCGGCGGCGCAG
GGGGGCGCCTGGCCCTGATCTGGATAGTCCCGCTCACCCTCAGCGGCCTCCTAGGAGTGGCCTGGGGGGC
ATCCAGTTTGGGAGCGCACCACATCCACCATTTCCATGGCAGCAGCAAGCATCATTCAGTGCCTATTGCA
ATCTACAGGTCACCAGCATCCTTGCGAGGCGGACACGCTGGGACAACATATATCTTTAGCAAAGGTGGTG
GACAGATTACATATAAGTGGCCTCCCAATGACCGCCCCAGTACACGAGCAGACAGGCTGGCCATCGGATT
TAGCACTGTCCAGAAGGAAGCCGTGTTGGTGCGTGTGGACAGTTCCTCAGGACTGGGTGACTACCTGGAG
CTGCACATACACCAAGGAAAAATTGGAGTTAAGTTTAATGTTGGGACAGATGACATCGCCATCGAAGAGT
CTAATGCAATCATTAATGATGGGAAATACCATGTAGTACGTTTCACAAGGAGTGGTGGCAATGCCACGTT
ACAGGTGGACAGCTGGCCAGTTATCGAACGCTACCCTGCAGGGCGTCAGCTCACAATCTTCAATAGCCAA
GCAACCATAATAATTGGCGGGAAAGAGCAGGGCCAGCCCTTCCAGGGCCAGCTCTCTGGTCTTTACTACA
ATGGCTTGAAAGTTCTGAATATGGCAGCAGAGAACGATGCCAACATCGCCATAGTGGGGAATGTGAGGCT
GGTCGGTGAAGTGCCTTCCTCTATGACAACTGAGTCGACAGCCACTGCCATGCAGTCTGAGATGTCCACC
TCAATCATGGAGACCACCACAACCCTGGCTACCAGCACAGCTAGAAGAGGCAAGCCCCCCACAAAGGAAC
CCATCAGCCAGACCACAGATGACATCCTTGTGGCCTCGGCAGAGTGTCCCAGTGACGATGAGGACATTGA
CCCCTGTGAGCCGAGCTCAGGTGGGTTAGCCAACCCCACCAGAGTAGGTGGCCGTGAACCATACCCAGGC
TCTGCAGAGGTGATTCGGGAGTCCAGCAGTACCACTGGCATGGTGGTGGGGATTGTCGCAGCAGCCGCTC
TGTGCATCCTCATCCTCCTCTATGCCATGTACAAGTACAGGAACCGGGATGAAGGGTCGTACCACGTGGA
TGAGAGTCGAAACTACATCAGTAACTCAGCACAGTCCAATGGGGCTGTGGTCAAGGAGAAGCAACCCAGC
AGTGCTAAAAGCGCCAACAAAAACAAGAAAAACAAGGATAAGGAGTATTACGTCTGA
""")
tools.deconstruct_dna_sequence(nrx1b, "NRX1B_SS4IN", True)

back_translated_sequence = tools.back_translate_amino_acid_sequence("GNFDNERLAIARQRIPYRLGRVVDEWLLDK")

tools.deconstruct_dna_sequence(back_translated_sequence, "NRX1B_IN2", True)

tools.insert_sequence_into_peptide_position("NRX1B_SS4IN", "NRX1B_IN2", "NRX1B_SS4IP2", 201)

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_deconstructed_sequences["NRX1B_SS4IP2"], 47, 392, "NRX1B_SS4IP2_ExC")

ecoR1Seq = tools.create_seq_object_from_string("GAATTC")
tools.deconstruct_dna_sequence(ecoR1Seq, "EcoRI", False)

seqPeptideSeq = tools.create_seq_object_from_string("ATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCT")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SecretionSignal", True)

linker1Seq = tools.create_seq_object_from_string("GGTACC")
tools.deconstruct_dna_sequence(linker1Seq, "Linker1", True)

linker2Seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker2Seq, "Linker2", True)

apTagSeq = tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstruct_dna_sequence(apTagSeq, "APTag", True)

linker3eq = tools.create_seq_object_from_string("GGTAGCGGA")
tools.deconstruct_dna_sequence(linker3eq, "Linker3", True)

tevSeq = tools.create_seq_object_from_string("GAGAACCTATACTTCCAAGGA")
tools.deconstruct_dna_sequence(tevSeq, "TEV", True)

gstSeq = tools.create_seq_object_from_string("ATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTTTTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAATGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGATGGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAACATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTGGATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTTGATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAAACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGATGTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAAAAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCATGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAA")
tools.deconstruct_dna_sequence(gstSeq, "GST", True)

hisTagSeq = tools.create_seq_object_from_string("CACCACCATCATCACCAC")
tools.deconstruct_dna_sequence(hisTagSeq, "HIS", True)

stopCodonsSeq = tools.create_seq_object_from_string("TAGTAA")
tools.deconstruct_dna_sequence(stopCodonsSeq, "STOPS", True)

# pd1ExtracellularPeptideSequence = tools.createSeqObjectFromString('LEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGM')

tools.create_construct_from_deconstructed_sequences(['SecretionSignal', 'Linker1', 'NRX1B_SS4IP2_ExC', "Linker2", 'APTag', 'Linker3', 'TEV', 'GST', 'HIS', 'STOPS'], 'NRX1B_SS4IP2_ExC-APTag-TEV-GST-HIS')
tools.create_construct_from_deconstructed_sequences(['NheI', 'SecretionSignal', 'Linker1', 'NRX1B_SS4IP2_ExC', "Linker2", 'APTag', 'Linker3', 'TEV', 'GST', 'HIS', 'STOPS', 'EcoRI'], 'NheI-NRX1B_SS4IP2_ExC-APTag-TEV-GST-HIS-EcoRI')

# pd1PeptideSequence = tools.createSeqObjectFromString('MWVRQVPWSFTWAVLQLSWQSGWLLEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGMVIGIMSALVGIPVLLLLAWALAVFCSTSMSEARGAGSKDDTLKEEPSAAPVPSVAYEELDFQGREKTPELPTACVHTEYATIVFTEGLGASAMGRRGSADGLQGPRPPRHEDGHCSWPL')


# compare1 = tools.comparePeptideConstructToSequence(tools.allConstructs['PD1_Extracellular'], pd1ExtracellularPeptideSequence)
# compare2 = tools.comparePeptideConstructToSequence(tools.allConstructs['NM_008798.2'], pd1PeptideSequence)

pass