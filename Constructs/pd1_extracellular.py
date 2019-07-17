from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.importSeqeuenceByNCBIIdentifier("NM_008798.2")
tools.deconstructImportedCDNASequence(tools.allSequences["NM_008798.2"], "NM_008798.2", 288)
tools.makeNewDeconstructedSequenceFromDeconstructedSequencePeptideRange(tools.allDeconstructedSequences["NM_008798.2"], 25, 169, "PD1_Extracellular")
dnaSeq = tools.returnDnaSequenceFromDeconstructedList(tools.allDeconstructedSequences["PD1_Extracellular"]['deconstructedList'])

nhe1Seq = tools.createSeqObjectFromString("GCTAGC")
tools.deconstructDNASequence(nhe1Seq, "NheI", False)

ecoR1Seq = tools.createSeqObjectFromString("GAATTC")
tools.deconstructDNASequence(ecoR1Seq, "EcoRI", False)

seqPeptideSeq = tools.createSeqObjectFromString("ATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCT")
tools.deconstructDNASequence(seqPeptideSeq, "SecretionSignal", True)

linker1Seq = tools.createSeqObjectFromString("GGTACC")
tools.deconstructDNASequence(linker1Seq, "Linker1", True)

linker2Seq = tools.createSeqObjectFromString("GGTAGTGGTGGTAGTGGT")
tools.deconstructDNASequence(linker2Seq, "Linker2", True)

apTagSeq = tools.createSeqObjectFromString("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")
tools.deconstructDNASequence(apTagSeq, "APTag", True)

linker3eq = tools.createSeqObjectFromString("GGTAGCGGA")
tools.deconstructDNASequence(linker3eq, "Linker3", True)

tevSeq = tools.createSeqObjectFromString("GAGAACCTATACTTCCAAGGA")
tools.deconstructDNASequence(tevSeq, "TEV", True)

hisTagSeq = tools.createSeqObjectFromString("CACCACCATCATCACCAC")
tools.deconstructDNASequence(hisTagSeq, "HIS", True)

stopCodonsSeq = tools.createSeqObjectFromString("TAGTAA")
tools.deconstructDNASequence(stopCodonsSeq, "STOPS", True)

pd1ExtracellularPeptideSequence = tools.createSeqObjectFromString('LEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGM')

tools.createConstructFromDeconstructedSequences(['SecretionSignal', 'Linker1', 'PD1_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS','STOPS'], 'PD1_Extracellular-APTag-TEV-HIS')
tools.createConstructFromDeconstructedSequences(['NheI' ,'SecretionSignal', 'Linker1', 'PD1_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS','STOPS', 'EcoRI'], 'NheI-PD1_Extracellular-APTag-TEV-HIS-EcoRI')

pd1PeptideSequence = tools.createSeqObjectFromString('MWVRQVPWSFTWAVLQLSWQSGWLLEVPNGPWRSLTFYPAWLTVSEGANATFTCSLSNWSEDLMLNWNRLSPSNQTEKQAAFCNGLSQPVQDARFQIIQLPNRHDFHMNILDTRRNDSGIYLCGAISLHPKAKIEESPGAELVVTERILETSTRYPSPSPKPEGRFQGMVIGIMSALVGIPVLLLLAWALAVFCSTSMSEARGAGSKDDTLKEEPSAAPVPSVAYEELDFQGREKTPELPTACVHTEYATIVFTEGLGASAMGRRGSADGLQGPRPPRHEDGHCSWPL')


compare1 = tools.comparePeptideConstructToSequence(tools.allConstructs['PD1_Extracellular'], pd1ExtracellularPeptideSequence)
compare2 = tools.comparePeptideConstructToSequence(tools.allConstructs['NM_008798.2'], pd1PeptideSequence)

pass