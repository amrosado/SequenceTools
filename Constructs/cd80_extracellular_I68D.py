from SequenceTools import SequenceTools

tools = SequenceTools(email="arosado@gatech.edu")
tools.importSeqeuenceByNCBIIdentifier("AJ278965.1")
tools.deconstructImportedCDNASequence(tools.allSequences["AJ278965.1"], "AJ278965.1", maxPeptideLength=306)
tools.makeNewDeconstructedSequenceFromConstructSequenceWithPeptideMutation(tools.allConstructs["AJ278965.1"], 'CD80_I68D', 68, 'I', 'D')
tools.makeNewDeconstructedSequenceFromDeconstructedSequencePeptideRange(tools.allConstructs["CD80_I68D"], 38, 246, "CD80_ExtracellularMTI68D")
dnaSeq = tools.returnDnaSequenceFromDeconstructedList(tools.allDeconstructedSequences["CD80_ExtracellularMTI68D"]['deconstructedList'])
dnaSeq2 = tools.createSeqObjectFromString("GTTTCCGTGGAGACGCAAGCTTATTTCAATGGGACTGCATATCTGCCGTGCCCAT TTACAAAGGCTCAAAACATAAGCCTGAGTGAGCTGGTAGTATTTTGGCAGGACCAGCAAA AGTTGGTTCTGTACGAGCACTATTTGGGCACAGAGAAACTTGATAGTGTGAATGCCAAGT ACCTGGGCCGCACGAGCTTTGACAGGAACAACTGGACTCTACGACTTCACAATGTTCAGA TCAAGGACATGGGCTCGTATGATTGTTTTATACAAAAAAAGCCACCCACAGGATCAATTA TCCTCCAACAGACATTAACAGAACTGTCAGTGATCGCCAACTTCAGTGAACCTGAAATAA AACTGGCTCAGAATGTAACAGGAAATTCTGGCATAAATTTGACCTGCACGTCTAAGCAAG GTCACCCGAAACCTAAGAAGATGTATTTTCTGATAACTAATTCAACTAATGAGTATGGTG ATAACATGCAGATATCACAAGATAATGTCACAGAACTGTTCAGTATCTCCAACAGCCTCT CTCTTTCATTCCCGGATGGTGTGTGGCATATGACCGTTGTGTGTGTTCTGGAAACGGAGT CAATGAAGATTTCCTCCAAACCTCTCAATTTCACTCAAGAGTTTCCATCTCCTCAAACGT ATTGGAAG")

nhe1Seq = tools.createSeqObjectFromString("GCTAGC")
tools.deconstructDNASequence(nhe1Seq, "NheI", False)

ecoR1Seq = tools.createSeqObjectFromString("GAATTC")
tools.deconstructDNASequence(ecoR1Seq, "EcoRI", False)

seqPeptideSeq = tools.createSeqObjectFromString("ATGGGGATCCTTCCCAGCCCTGGGATGCCTGCGCTGCTCTCCCTCGTGAGCCTTCTCTCCGTGCTGCTGATGGGTTGCGTAGCT")
tools.deconstructDNASequence(seqPeptideSeq, "SecretionSignal", True)

linker1Seq = tools.createSeqObjectFromString("GGTACCGGA")
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

cd80I68DExtracellularPeptideSequence = tools.createSeqObjectFromString('VDEQLSKSVKDKVLLPCRYNSPHEDESEDRDYWQKHDKVVLSVIAGKLKVWPEYKNRTLYDNTTYSLIILGLVLSDRGTYSCVVQKKERGTYEVKHLALVKLSIKADFSTPNITESGNPSADTKRITCFASGGFPKPRFSWLENGRELPGINTTISQDPESELYTISSQLDFNTTRNHTIKCLIKYGDAHVSEDFTWEKPPEDPPDSKN')

tools.createConstructFromDeconstructedSequences(['SecretionSignal', 'Linker1', 'CD80_ExtracellularMTI68D', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS','STOPS'], 'CD80_ExtracellularMTI68D-APTag-TEV-HIS')
tools.createConstructFromDeconstructedSequences(['NheI' ,'SecretionSignal', 'Linker1', 'CD80_ExtracellularMTI68D', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS','STOPS', 'EcoRI'], 'NheI-CD80_ExtracellularMTI68D-APTag-TEV-HIS-EcoRI')

compare1 = tools.comparePeptideConstructToSequence(tools.allConstructs['CD80_ExtracellularMTI68D'], cd80I68DExtracellularPeptideSequence)

pass