from SequenceTools import SequenceTools

tools = SequenceTools()
tools.importSeqeuenceByNCBIIdentifier("NM_009843.4")
tools.deconstructImportedCDNASequence(tools.allSequences["NM_009843.4"], "NM_009843.4", minPeptideLength=223)

pass