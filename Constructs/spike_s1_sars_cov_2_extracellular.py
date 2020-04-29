from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.import_sequence_by_ncbi_identifier("NC_045512.2")
tools.deconstruct_imported_orf_sequence(tools.all_sequences["NC_045512.2"], "NC_045512.2", 'MFVFLVLLPL',  min_peptide_length=1273)
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_deconstructed_sequences["NC_045512.2"], 13, 685, "SarsCov2SpikeS1")

notISeq = tools.create_seq_object_from_string("GCGGCCGC")
tools.deconstruct_dna_sequence(notISeq, "NotI", False)

xbaISeq = tools.create_seq_object_from_string("TCTAGA")
tools.deconstruct_dna_sequence(xbaISeq, "XbaI", False)

seqPeptideSeq = tools.create_seq_object_from_string("ATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCT")
tools.deconstruct_dna_sequence(seqPeptideSeq, "SecretionSignal", True)

linker1Seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker1Seq, "Linker1", True)

linker2Seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker2Seq, "Linker2", True)

linker3eq = tools.create_seq_object_from_string("GGTAGCGGA")
tools.deconstruct_dna_sequence(linker3eq, "Linker3", True)

apTagSeq = tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstruct_dna_sequence(apTagSeq, "APTag", True)

tevSeq = tools.create_seq_object_from_string("GAGAACCTATACTTCCAAGGA")
tools.deconstruct_dna_sequence(tevSeq, "TEV", True)

hisTagSeq = tools.create_seq_object_from_string("CACCACCATCATCACCAC")
tools.deconstruct_dna_sequence(hisTagSeq, "HIS", True)

stopCodonsSeq = tools.create_seq_object_from_string("TAGTAA")
tools.deconstruct_dna_sequence(stopCodonsSeq, "STOPS", True)

tools.create_construct_from_deconstructed_sequences(['SecretionSignal', 'Linker1', 'SarsCov2SpikeS1', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS'], 'Sec_SarsCov2SpikeS1_BAP_TEV_HIS')

tools.create_construct_from_deconstructed_sequences(['NotI', 'SecretionSignal', 'Linker1', 'SarsCov2SpikeS1', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS', 'STOPS', 'XbaI'], 'NotI-SarsCov2SpikeS1-APTag-TEV-HIS-XbaI')

sars_cov_2_spike = tools.create_seq_object_from_string(
    '''
    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS
NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV
NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE
GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT
LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK
CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN
CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD
YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC
NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN
FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP
GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY
ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI
SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE
VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC
LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM
QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN
TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA
SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA
ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP
LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL
QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD
SEPVLKGVKLHYT
    '''
)

sars_cov_2_spike_s1 = tools.create_seq_object_from_string(
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

compare_1 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['SarsCov2SpikeS1'], sars_cov_2_spike_s1)
compare_2 = tools.compare_peptide_construct_to_sequence(tools.all_constructs['NC_045512.2'], sars_cov_2_spike)
pass