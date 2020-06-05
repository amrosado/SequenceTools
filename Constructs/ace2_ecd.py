from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("NM_001371415.1")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NM_001371415.1"], "NM_001371415.1", 'MSSSSWLLLS',  min_peptide_length=492)

tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_constructs["NM_001371415.1"], 18, 740, "ACE2_ECD")

notISeq = tools.create_seq_object_from_string("GCGGCCGC")

tools.deconstruct_dna_sequence(notISeq, "NotI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")

tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

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

seqPeptideSeq = tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SecretionSignal_mouseIgKappa", True)


ha3xpeptide = tools.create_seq_object_from_string(
    """
    MEYPYDVPDYAAEYPYDVPDYAAEYPYDVPDYAAKLE
    """
)

ha3x_dna_seq = tools.back_translate_amino_acid_sequence_random_codon(ha3xpeptide)

tools.deconstruct_dna_sequence(ha3x_dna_seq, "Ha3x", True)

ha3x = tools.create_seq_object_from_string(
    """
    tacccatacgatgttccagattacgcttatccttatgacgtacctgactatgcatacccttatgatgtaccagactacgct
    """
)

tools.deconstruct_dna_sequence(ha3x, "3xHA", True)

flag3x = tools.create_seq_object_from_string(
    '''
    gactacaaagaccatgacggtgattataaagatcatgacatcgattacaaggatgacgatgacaag
    '''
)

tools.deconstruct_dna_sequence(flag3x, 'FLAG3x', True)

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

linker_3_mclellan = tools.create_seq_object_from_string(
    '''
    GGC
    '''
)

tools.deconstruct_dna_sequence(linker_3_mclellan, 'Linker_3', True)


construct_list_without_re = ['SecretionSignal_mouseIgKappa', 'ACE2_ECD', 'Linker_2', 'APTag', 'Linker_0', 'Ha3x', "Linker_1", 'hrv3c_protease_cleavage', 'Linker_4', 'his_8', 'stops']

construct_list_with_re = construct_list_without_re.copy()

construct_list_with_re.append('XbaI')
construct_list_with_re.insert(0, 'NotI')

tools.create_construct_from_deconstructed_sequences(construct_list_without_re, 'ACE2_Ha3x_BAP_HRV3C_His8x')
tools.create_construct_from_deconstructed_sequences(construct_list_with_re, 'NotI_ACE2_Ha3x_BAP_HRV3C_His8x_XbaI')

ace2_ecd_uniprot = tools.create_seq_object_from_string(
    '''
    QSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQS
TLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDN
PQECLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYE
DYGDYWRGDYEVNGVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYIS
PIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVS
VGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMG
HIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEIN
FLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETY
CDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNM
LRLGKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQS
IKVRISLKSALGDKAYEWNDNEMYLFRSSVAYAMRQYFLKVKNQMILFGEEDVRVANLKP
RISFNFFVTAPKNVSDIIPRTEVEKAIRMSRSRINDAFRLNDNSLEFLGIQPTLGPPNQP
PVS
    '''
)

compare_3xha = tools.compare_peptide_construct_to_sequence(tools.all_constructs['Ha3x'], ha3xpeptide)

compare_ace2_ecd = tools.compare_peptide_construct_to_sequence(tools.all_constructs['ACE2_ECD'], ace2_ecd_uniprot)

working_seq = tools.create_seq_object_from_string(
    """
    atggagacagacacactcctgctatgggtactgctgctctgggttccaggttccactggtg
    accagtccaccattgaggaacaggccaagacatttttggacaagtttaaccacgaagccga
    agacctgttctatcaaagttcacttgcttcttggaattataacaccaatattactgaagag
    aatgtccaaaacatgaataatgctggggacaaatggtctgcctttttaaaggaacagtcca
    cacttgcccaaatgtatccactacaagaaattcagaatctcacagtcaagcttcagctgca
    ggctcttcagcaaaatgggtcttcagtgctctcagaagacaagagcaaacggttgaacaca
    attctaaatacaatgagcaccatctacagtactggaaaagtttgtaacccagataatccac
    aagaatgcttattacttgaaccaggtttgaatgaaataatggcaaacagtttagactacaa
    tgagaggctctgggcttgggaaagctggagatctgaggtcggcaagcagctgaggccatta
    tatgaagagtatgtggtcttgaaaaatgagatggcaagagcaaatcattatgaggactatg
    gggattattggagaggagactatgaagtaaatggggtagatggctatgactacagccgcgg
    ccagttgattgaagatgtggaacatacctttgaagagattaaaccattatatgaacatctt
    catgcctatgtgagggcaaagttgatgaatgcctatccttcctatatcagtccaattggat
    gcctccctgctcatttgcttggtgatatgtggggtagattttggacaaatctgtactcttt
    gacagttccctttggacagaaaccaaacatagatgttactgatgcaatggtggaccaggcc
    tgggatgcacagagaatattcaaggaggccgagaagttctttgtatctgttggtcttccta
    atatgactcaaggattctgggaaaattccatgctaacggacccaggaaatgttcagaaagc
    agtctgccatcccacagcttgggacctggggaagggcgacttcaggatccttatgtgcaca
    aaggtgacaatggacgacttcctgacagctcatcatgagatggggcatatccagtatgata
    tggcatatgctgcacaaccttttctgctaagaaatggagctaatgaaggattccatgaagc
    tgttggggaaatcatgtcactttctgcagccacacctaagcatttaaaatccattggtctt
    ctgtcacccgattttcaagaagacaatgaaacagaaataaacttcctgctcaaacaagcac
    tcacgattgttgggactctgccatttacttacatgttagagaagtggaggtggatggtctt
    taaaggggaaattcccaaagaccagtggatgaaaaagtggtgggagatgaagcgagagata
    gttggggtggtggaacctgtgccccatgatgaaacatactgtgaccccgcatctctgttcc
    atgtttctaatgattactcattcattcgatattacacaaggaccctttaccaattccagtt
    tcaagaagcactttgtcaagcagctaaacatgaaggccctctgcacaaatgtgacatctca
    aactctacagaagctggacagaaactgttcaatatgctgaggcttggaaaatcagaaccct
    ggaccctagcattggaaaatgttgtaggagcaaagaacatgaatgtaaggccactgctcaa
    ctactttgagcccttatttacctggctgaaagaccagaacaagaattcttttgtgggatgg
    agtaccgactggagtccatatgcagaccaaagcatcaaagtgaggataagcctaaaatcag
    ctcttggagataaagcatatgaatggaacgacaatgaaatgtacctgttccgatcatctgt
    tgcatatgctatgaggcagtactttttaaaagtaaaaaatcagatgattctttttggggag
    gaggatgtgcgagtggctaatttgaaaccaagaatctcctttaatttctttgtcactgcac
    ctaaaaatgtgtctgatatcattcctagaactgaagttgaaaaggccatcaggatgtcccg
    gagccgtatcaatgatgctttccgtctgaatgacaacagcctagagtttctggggatacag
    ccaacacttggacctcctaaccagccccctgtttccagatccggtctgaatgatattttcg
    aagcgcagaaaattgaatggcatgaaggtaccggaatggaatacccatatgatgtaccaga
    ctacgctgcagaatatccctacgatgttccagactatgccgcggagtacccctatgatgta
    cccgattacgccgcaaaattggaaggatccggactggaggtgctgttccagggcccaagcg
    cccatcaccaccatcaccaccatcattgataatga
    """
)

pass