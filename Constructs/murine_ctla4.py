from SequenceTools import SequenceTools

tools = SequenceTools()
tools.import_sequence_by_ncbi_identifier("NM_009843.4")
tools.deconstruct_imported_cdna_sequence(tools.all_sequences["NM_009843.4"], "NM_009843.4", min_peptide_length=223)

pass