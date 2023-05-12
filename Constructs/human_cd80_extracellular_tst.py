from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("AJ278965.1")
tools.deconstruct_imported_cdna_sequence(tools.all_sequences["AJ278965.1"], "AJ278965.1", maxPeptideLength=306)
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_deconstructed_sequences["AJ278965.1"], 38, 246, "CD80_Extracellular")

dnaSeq = tools.return_dna_sequence_from_deconstructed_list(tools.all_deconstructed_sequences["CD80_Extracellular"]['deconstructedList'])
dnaSeq2 = tools.create_seq_object_from_string("GTTTCCGTGGAGACGCAAGCTTATTTCAATGGGACTGCATATCTGCCGTGCCCAT TTACAAAGGCTCAAAACATAAGCCTGAGTGAGCTGGTAGTATTTTGGCAGGACCAGCAAA AGTTGGTTCTGTACGAGCACTATTTGGGCACAGAGAAACTTGATAGTGTGAATGCCAAGT ACCTGGGCCGCACGAGCTTTGACAGGAACAACTGGACTCTACGACTTCACAATGTTCAGA TCAAGGACATGGGCTCGTATGATTGTTTTATACAAAAAAAGCCACCCACAGGATCAATTA TCCTCCAACAGACATTAACAGAACTGTCAGTGATCGCCAACTTCAGTGAACCTGAAATAA AACTGGCTCAGAATGTAACAGGAAATTCTGGCATAAATTTGACCTGCACGTCTAAGCAAG GTCACCCGAAACCTAAGAAGATGTATTTTCTGATAACTAATTCAACTAATGAGTATGGTG ATAACATGCAGATATCACAAGATAATGTCACAGAACTGTTCAGTATCTCCAACAGCCTCT CTCTTTCATTCCCGGATGGTGTGTGGCATATGACCGTTGTGTGTGTTCTGGAAACGGAGT CAATGAAGATTTCCTCCAAACCTCTCAATTTCACTCAAGAGTTTCCATCTCCTCAAACGT ATTGGAAG")

xba_1 = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xba_1, "XbaI", False)

xho_1 = tools.create_seq_object_from_string("CTCGAG")
tools.deconstruct_dna_sequence(xho_1, "XhoI", False)

kozak_seq_non_coding = tools.create_seq_object_from_string("gccAcc")
tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakNonCoding", False)

kozak_seq_coding = tools.create_seq_object_from_string("ATGGGT")
tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakCoding", True)

# Secretion signal mouse Ig kappa chain
seqPeptideSeq = tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SS_MIgK", True)

linker1Seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker1Seq, "Linker1", True)

linker2Seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker2Seq, "Linker2", True)

apTagSeq = tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstruct_dna_sequence(apTagSeq, "APTag", True)

linker3eq = tools.create_seq_object_from_string("GGTAGCGGA")
tools.deconstruct_dna_sequence(linker3eq, "Linker3", True)

tevSeq = tools.create_seq_object_from_string("GAGAACCTATACTTCCAAGGA")
tools.deconstruct_dna_sequence(tevSeq, "TEV", True)

twin_strept_2_tag = tools.create_seq_object_from_string(
    '''
    TGGTCCCACCCCCAGTTCGAGAAGGGCGGCGGTAGTGGAGGGGGCGGATCTGGCGGCTCAGCTTGGAGCCACCCCCAGTTCGAAAAG
    '''
)
tools.deconstruct_dna_sequence(twin_strept_2_tag, "TST", True)

linker_4 = tools.create_seq_object_from_string(
    '''
    GGC
    '''
)
tools.deconstruct_dna_sequence(linker_4, "Linker4", True)

linker_5 = tools.create_seq_object_from_string(
    '''
    AGCGCC
    '''
)
tools.deconstruct_dna_sequence(linker_5, "TST", True)

eight_x_his_tag = tools.create_seq_object_from_string(
    '''
    CATCACCACCATCACCACCATCAT
    '''
)
tools.deconstruct_dna_sequence(eight_x_his_tag, "8xHis", True)

ires2_egfp = tools.create_seq_object_from_string(
    '''
    cccctctccctcccccccccctaacgttactggccgaagccgcttggaataaggccggtgtgcgtttgtctatatgttattttccaccatattgcc
    gtcttttggcaatgtgagggcccggaaacctggccctgtcttcttgacgagcattcctaggggtctttcccctctcgccaaaggaatgcaaggtct
    gttgaatgtcgtgaaggaagcagttcctctggaagcttcttgaagacaaacaacgtctgtagcgaccctttgcaggcagcggaaccccccacctgg
    cgacaggtgcctctgcggccaaaagccacgtgtataagatacacctgcaaaggcggcacaaccccagtgccacgttgtgagttggatagttgtgga
    aagagtcaaatggctctcctcaagcgtattcaacaaggggctgaaggatgcccagaaggtaacccattgtatgggatctgatctggggcctcggta
    cacatgctttacatgtgtttagtcgaggttaaaaaaacgtctaggccccccgaaccacggggacgtggttttcctttgaaaaacacgatgataata
    tggccacaaccatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttca
    gcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccc
    tcgtgaccaccttgacctacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggct
    acgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcg
    agctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaaggtctatatcaccgccgaca
    agcagaagaacggcatcaaggtgaacttcaagacccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacaccccca
    tcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcc
    tgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa
    '''
)

tools.deconstruct_dna_sequence(ires2_egfp, "IRES2_EGFP", True)

stopCodonsSeq = tools.create_seq_object_from_string("TAGTAA")
tools.deconstruct_dna_sequence(stopCodonsSeq, "STOPS", True)

cd80ExtracellularPeptideSequence = tools.create_seq_object_from_string('VDEQLSKSVKDKVLLPCRYNSPHEDESEDRIYWQKHDKVVLSVIAGKLKVWPEYKNRTLYDNTTYSLIILGLVLSDRGTYSCVVQKKERGTYEVKHLALVKLSIKADFSTPNITESGNPSADTKRITCFASGGFPKPRFSWLENGRELPGINTTISQDPESELYTISSQLDFNTTRNHTIKCLIKYGDAHVSEDFTWEKPPEDPPDSKN')

tools.create_construct_from_deconstructed_sequences(['KozakNonCoding', 'KozakCoding', 'SS_MIgK', 'Linker1', 'CD80_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'TST', 'Linker4', '8xHis'], 'CD80Extracellular-APTag-TEV-TST-8xHis-IRES_EGFP')
tools.create_construct_from_deconstructed_sequences(['XbaI' , 'KozakNonCoding', 'KozakCoding', 'SS_MIgK', 'Linker1', 'CD80_Extracellular', "Linker2", 'APTag', 'Linker3', 'TST', 'Linker4', '8xHis', 'XhoI'], 'XbaI-CD80Extracellular-APTag-TEV-TST-8xHis-IRES_EGFP-XhoI')

cd80PeptideSequence = tools.create_seq_object_from_string('MACNCQLMQDTPLLKFPCPRLILLFVLLIRLSQVSSDVDEQLSKSVKDKVLLPCRYNSPHEDESEDRIYWQKHDKVVLSVIAGKLKVWPEYKNRTLYDNTTYSLIILGLVLSDRGTYSCVVQKKERGTYEVKHLALVKLSIKADFSTPNITESGNPSADTKRITCFASGGFPKPRFSWLENGRELPGINTTISQDPESELYTISSQLDFNTTRNHTIKCLIKYGDAHVSEDFTWEKPPEDPPDSKNTLVLFGAGFGAVITVVVIVVIIKCFCKHRSCFRRNEASRETNNSLTFGPEEALAEQTVFL')
genScriptSeq = tools.create_seq_object_from_string('''
GCTAGCATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCTGGTAC
CGGAGTTGATGAACAACTGTCCAAGTCAGTGAAAGATAAGGTATTGCTGCCTTGCCGTTACAACTCTCCTCATGAAGATGAGTCTGAAGACCGAATCTAC
TGGCAAAAACATGACAAAGTGGTGCTGTCTGTCATTGCTGGGAAACTAAAAGTGTGGCCCGAGTATAAGAACCGGACTTTATATGACAACACTACCTACT
CTCTTATCATCCTGGGCCTGGTCCTTTCAGACCGGGGCACATACAGCTGTGTCGTTCAAAAGAAGGAAAGAGGAACGTATGAAGTTAAACACTTGGCTTT
AGTAAAGTTGTCCATCAAAGCTGACTTCTCTACCCCCAACATAACTGAGTCTGGAAACCCATCTGCAGACACTAAAAGGATTACCTGCTTTGCTTCCGGG
GGTTTCCCAAAGCCTCGCTTCTCTTGGTTGGAAAATGGAAGAGAATTACCTGGCATCAATACGACAATTTCCCAGGATCCTGAATCTGAATTGTACACCA
TTAGTAGCCAACTAGATTTCAATACGACTCGCAACCACACCATTAAGTGTCTCATTAAATATGGAGATGCTCACGTGTCAGAGGACTTCACCTGGGAAAA
ACCCCCAGAAGACCCTCCTGATAGCAAGAACGGTAGTGGTGGTAGTGGTGGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAAGGTAGC
GGAGAGAACCTATACTTCCAAGGACACCACCATCATCACCACTAGTAAGAATTC
''')

genScriptSeq2 = tools.create_seq_object_from_string('''
GCTAGCATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCTGGTAC
CGGAGTTGATGAACAACTGTCCAAGTCAGTGAAAGATAAGGTATTGCTGCCTTGCCGTTACAACTCTCCTCATGAAGATGAGTCTGAAGACCGAATCTAC
TGGCAAAAACATGACAAAGTGGTGCTGTCTGTCATTGCTGGGAAACTAAAAGTGTGGCCCGAGTATAAGAACCGGACTTTATATGACAACACTACCTACT
CTCTTATCATCCTGGGCCTGGTCCTTTCAGACCGGGGCACATACAGCTGTGTCGTTCAAAAGAAGGAAAGAGGAACGTATGAAGTTAAACACTTGGCTTT
AGTAAAGTTGTCCATCAAAGCTGACTTCTCTACCCCCAACATAACTGAGTCTGGAAACCCATCTGCAGACACTAAAAGGATTACCTGCTTTGCTTCCGGG
GGTTTCCCAAAGCCTCGCTTCTCTTGGTTGGAAAATGGAAGAGAATTACCTGGCATCAATACGACAATTTCCCAGGATCCTGAATCTGAATTGTACACCA
TTAGTAGCCAACTAGATTTCAATACGACTCGCAACCACACCATTAAGTGTCTCATTAAATATGGAGATGCTCACGTGTCAGAGGACTTCACCTGGGAAAA
ACCCCCAGAAGACCCTCCTGATAGCAAGAACGGTAGTGGTGGTAGTGGTGGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAAGGTAGC
GGAGAGAACCTATACTTCCAAGGACACCACCATCATCACCACTAGTAAGAATTC
''')

genScriptSeq3 = tools.create_seq_object_from_string('''
GCTAGCATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCTGGTAC
CGGAGTTGATGAACAACTGTCCAAGTCAGTGAAAGATAAGGTATTGCTGCCTTGCCGTTACAACTCTCCTCATGAAGATGAGTCTGAAGACCGAGACTAC
TGGCAAAAACATGACAAAGTGGTGCTGTCTGTCATTGCTGGGAAACTAAAAGTGTGGCCCGAGTATAAGAACCGGACTTTATATGACAACACTACCTACT
CTCTTATCATCCTGGGCCTGGTCCTTTCAGACCGGGGCACATACAGCTGTGTCGTTCAAAAGAAGGAAAGAGGAACGTATGAAGTTAAACACTTGGCTTT
AGTAAAGTTGTCCATCAAAGCTGACTTCTCTACCCCCAACATAACTGAGTCTGGAAACCCATCTGCAGACACTAAAAGGATTACCTGCTTTGCTTCCGGG
GGTTTCCCAAAGCCTCGCTTCTCTTGGTTGGAAAATGGAAGAGAATTACCTGGCATCAATACGACAATTTCCCAGGATCCTGAATCTGAATTGTACACCA
TTAGTAGCCAACTAGATTTCAATACGACTCGCAACCACACCATTAAGTGTCTCATTAAATATGGAGATGCTCACGTGTCAGAGGACTTCACCTGGGAAAA
ACCCCCAGAAGACCCTCCTGATAGCAAGAACGGTAGTGGTGGTAGTGGTGGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAAGGTAGC
GGAGAGAACCTATACTTCCAAGGACACCACCATCATCACCACTAGTAAGAATTC
''')

# compare1 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['CD80_Extracellular'], cd80ExtracellularPeptideSequence)
# compare2 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['AJ278965.1'], cd80PeptideSequence)
# compareGenScript = tools.compare_dna_construct_to_sequence(tools.all_constructs['NheI-CD80Extracellular-APTag-TEV-HIS-EcoRI'], genScriptSeq)
# compareGenScript2 = tools.compare_dna_construct_to_sequence(tools.all_constructs['NheI-CD80Extracellular-APTag-TEV-HIS-EcoRI'], genScriptSeq2)
# compareGenScript3 = tools.compare_dna_construct_to_sequence(tools.all_constructs['NheI-CD80Extracellular-APTag-TEV-HIS-EcoRI'], genScriptSeq3)

# translated = genScriptSeq2.translate()

pass