from SequenceTools import SequenceTools

tools = SequenceTools()

bira_seq = tools.create_seq_object_from_string(
    """
    atgaaggataacaccgtgccactgaaattgattgccctgttagcgaacggtgaatttcactctggcgagcagttgggtgaaacgctgggaatgagccgggcggctattaataaacacattcagacactgcgtgactggggcgttgatgtctttaccgttccgggtaaaggatacagcctgcctgagcctatccagttacttaatgctaaacagatattgggtcagctggatggcggtagtgtagccgtgctgcctgtgattgactccacgaatcagtaccttcttgatcgtatcggagagcttaaatcgggcgatgcttgcattgcagaataccagcaggctggccgtggtcgccggggtcggaaatggttttcgccttttggcgcaaacttatatttgtcgatgttctggcgtctggaacaaggcccggcggcggcgattggtttaagtctggttatcggtatcgtgatggcggaagtattacgcaagctgggtgcagataaagttcgtgttaaatggcctaatgacctctatctgcaggatcgcaagctggcaggcattctggtggagctgactggcaaaactggcgatgcggcgcaaatagtcattggagccgggatcaacatggcaatgcgccgtgttgaagagagtgtcgttaatcaggggtggatcacgctgcaggaagcggggatcaatctcgatcgtaatacgttggcggccatgctaatacgtgaattacgtgctgcgttggaactcttcgaacaagaaggattggcaccttatctgtcgcgctgggaaaagctggataattttattaatcgcccagtgaaacttatcattggtgataaagaaatatttggcatttcacgcggaatagacaaacagggggctttattacttgagcaggatggaataataaaaccctggatgggcggtgaaatatccctgcgtagtgcagaaaaa
    """
)

tools.deconstruct_dna_sequence(bira_seq, "BirA", True)

gst_seq = tools.create_seq_object_from_string(
    """
    ATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTTTTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAATGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGATGGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAACATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTGGATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTTGATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAAACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGATGTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAAAAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCATGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAA
    """
)

tools.deconstruct_dna_sequence(gst_seq, "GST", True)

gst_bira = tools.create_seq_object_from_string(
"""
atggctagcatgactggtggacagcaaatgggtcgcggatccgaattcATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTTTTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAATGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGATGGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAACATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTGGATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTTGATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAAACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGATGTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAAAAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCATGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAATCGGATCTGGTTCCGCGTGGATCCatgaaggataacaccgtgccactgaaattgattgccctgttagcgaacggtgaatttcactctggcgagcagttgggtgaaacgctgggaatgagccgggcggctattaataaacacattcagacactgcgtgactggggcgttgatgtctttaccgttccgggtaaaggatacagcctgcctgagcctatccagttacttaatgctaaacagatattgggtcagctggatggcggtagtgtagccgtgctgcctgtgattgactccacgaatcagtaccttcttgatcgtatcggagagcttaaatcgggcgatgcttgcattgcagaataccagcaggctggccgtggtcgccggggtcggaaatggttttcgccttttggcgcaaacttatatttgtcgatgttctggcgtctggaacaaggcccggcggcggcgattggtttaagtctggttatcggtatcgtgatggcggaagtattacgcaagctgggtgcagataaagttcgtgttaaatggcctaatgacctctatctgcaggatcgcaagctggcaggcattctggtggagctgactggcaaaactggcgatgcggcgcaaatagtcattggagccgggatcaacatggcaatgcgccgtgttgaagagagtgtcgttaatcaggggtggatcacgctgcaggaagcggggatcaatctcgatcgtaatacgttggcggccatgctaatacgtgaattacgtgctgcgttggaactcttcgaacaagaaggattggcaccttatctgtcgcgctgggaaaagctggataattttattaatcgcccagtgaaacttatcattggtgataaagaaatatttggcatttcacgcggaatagacaaacagggggctttattacttgagcaggatggaataataaaaccctggatgggcggtgaaatatccctgcgtagtgcagaaaaa
"""
)

tools.deconstruct_dna_sequence(gst_bira, "T7leader_GST_BirA", True)

hind3 = tools.create_seq_object_from_string(
    """
    aagctt
    """
)

bamhI = tools.create_seq_object_from_string(
    """
    ggatcc
    """
)

sacI = tools.create_seq_object_from_string(
    """
    gagctc
    """
)

ecoRI = tools.create_seq_object_from_string(
    """
    gaattc
    """
)

stopCodonsSeq = tools.create_seq_object_from_string("TAGTAA")
tools.deconstruct_dna_sequence(stopCodonsSeq, "STOPS", True)

tools.deconstruct_dna_sequence(hind3, "HindIII", False)

tools.deconstruct_dna_sequence(ecoRI, "EcoRI", False)

tools.deconstruct_dna_sequence(bamhI, "BamHI", False)

tools.deconstruct_dna_sequence(sacI, "SacI", False)

linker1Seq = tools.create_seq_object_from_string("GGTACCGGA")

tools.deconstruct_dna_sequence(linker1Seq, "Linker1", True)

construct_list_without_re = ["T7leader_GST_BirA"]

construct_list_without_re_2 = ["GST", "Linker1", "BirA", "STOPS"]

construct = construct_list_without_re.copy()

construct.append("HindIII")
construct.insert(0, "EcoRI")

constuct2 = construct_list_without_re.copy()

constuct2.append("HindIII")
constuct2.insert(0, "SacI")

construct_3 = construct_list_without_re_2.copy()

construct_3.append("HindIII")
construct_3.insert(0, "SacI")

# tools.create_construct_from_deconstructed_sequences(construct, 'EcoRI_T7leader_GST_thrombin_BiotinLigaseBirA_HindIII')

tools.create_construct_from_deconstructed_sequences(constuct2, 'SacI_T7leader_GST_thrombin_BiotinLigaseBirA_HindIII')

tools.create_construct_from_deconstructed_sequences(construct_3, "SacI_GST_BirA_HindIII")

pass