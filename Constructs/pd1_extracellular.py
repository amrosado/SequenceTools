from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("NM_008798.2")
tools.deconstruct_imported_cdna_sequence(tools.allSequences["NM_008798.2"], "NM_008798.2", 288)
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.allDeconstructedSequences["NM_008798.2"], 25, 169, "PD1_Extracellular")
dnaSeq = tools.return_dna_sequence_from_deconstructed_list(tools.allDeconstructedSequences["PD1_Extracellular"]['deconstructedList'])

nhe1Seq = tools.create_seq_object_from_string("GCTAGC")
tools.deconstruct_dna_sequence(nhe1Seq, "NheI", False)

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

hisTagSeq = tools.create_seq_object_from_string("CACCACCATCATCACCAC")
tools.deconstruct_dna_sequence(hisTagSeq, "HIS", True)

stopCodonsSeq = tools.create_seq_object_from_string("TAGTAA")
tools.deconstruct_dna_sequence(stopCodonsSeq, "STOPS", True)

pd1ExtracellularPeptideSequence = tools.create_seq_object_from_string('LEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGM')

tools.create_construct_from_deconstructed_sequences(['SecretionSignal', 'Linker1', 'PD1_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS', 'STOPS'], 'PD1_Extracellular-APTag-TEV-HIS')
tools.create_construct_from_deconstructed_sequences(['NheI' , 'SecretionSignal', 'Linker1', 'PD1_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS', 'STOPS', 'EcoRI'], 'NheI-PD1_Extracellular-APTag-TEV-HIS-EcoRI')

pd1PeptideSequence = tools.create_seq_object_from_string('MWVRQVPWSFTWAVLQLSWQSGWLLEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGMVIGIMSALVGIPVLLLLAWALAVFCSTSMSEARGAGSKDDTLKEEPSAAPVPSVAYEELDFQGREKTPELPTACVHTEYATIVFTEGLGASAMGRRGSADGLQGPRPPRHEDGHCSWPL')


compare1 = tools.compare_peptide_construct_to_sequence(tools.allConstructs['PD1_Extracellular'], pd1ExtracellularPeptideSequence)
compare2 = tools.compare_peptide_construct_to_sequence(tools.allConstructs['NM_008798.2'], pd1PeptideSequence)

pass