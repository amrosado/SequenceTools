from SequenceTools import SequenceTools
from CommonSeqsTags import CommonSeqsTags

common_tags = CommonSeqsTags()

tools = SequenceTools(email="arosado@gatech.edu")

tools.import_sequence_by_ncbi_identifier("NM_005191")
tools.deconstruct_imported_cdna_sequence(tools.all_sequences["NM_005191"], "NM_005191", min_peptide_length=288)
tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(tools.all_deconstructed_sequences["NM_005191"], 35, 242, "hCD80_Extracellular")

tools.import_deconstructed_sequence(common_tags.return_xba1())
tools.import_deconstructed_sequence(common_tags.return_xho1())
tools.import_deconstructed_sequence(common_tags.return_kozak_non_coding())
tools.import_deconstructed_sequence(common_tags.return_tst_tag())
tools.import_deconstructed_sequence(common_tags.return_higg1_fc())
tools.import_deconstructed_sequence(common_tags.return_ap_tag())
tools.import_deconstructed_sequence(common_tags.return_tev())
tools.import_deconstructed_sequence(common_tags.return_8xhis_tag())
tools.import_deconstructed_sequence(common_tags.return_stops())
tools.import_deconstructed_sequence(common_tags.return_ssmigk())

linker1Seq = tools.create_seq_object_from_string("GGTACCGGA")
tools.deconstruct_dna_sequence(linker1Seq, "Linker1", True)

linker2Seq = tools.create_seq_object_from_string("GGTAGTGGTGGTAGTGGT")
tools.deconstruct_dna_sequence(linker2Seq, "Linker2", True)

linker3eq = tools.create_seq_object_from_string("GGTAGCGGA")
tools.deconstruct_dna_sequence(linker3eq, "Linker3", True)

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

construct_list = ['KozakNonCoding', 'SSMIGGK', 'Linker1', 'hCD80_Extracellular',  "hIgG1_Fc", "Linker2", 'APTag', 'Linker3', 'TEV', '8xHis', 'TST', 'Stops']

cuts_construct_list = ['XbaI'] + construct_list

tools.create_construct_from_deconstructed_sequences(construct_list, 'hCD80Extracellular-hIgG1Fc-APTag-TEV-8xHis-TST')
tools.create_construct_from_deconstructed_sequences(cuts_construct_list, 'XbaI-hCD80Extracellular-hIgG1Fc-APTag-TEV-8xHis-TST')

hcd80ec_seq = tools.create_seq_object_from_string("VIHVTKEVKEVATLSCGHNVSVEELAQTRIYWQKEKKMVLTMMSGDMNIWPEYKNRTIFDITNNLSIVILALRPSDEGTYECVVLKYEKDAFKREHLAEVTLSVKADFPTPSISDFEIPTSNIRRIICSTSGGFPEPHLSWLENGEELNAINTTVSQDPETELYAVSSKLDFNMTTNHSFMCLIKYGHLRVNQTFNWNTTKQEHFPDN")
compare_result = tools.compare_peptide_construct_to_sequence(tools.all_constructs["hCD80_Extracellular"], hcd80ec_seq)

print("Print constructed DNA sequence for XbaI-hCD80Extracellular-hIgG1Fc-APTag-TEV-8xHis-TST\n")
print("{}".format(tools.all_constructs["XbaI-hCD80Extracellular-hIgG1Fc-APTag-TEV-8xHis-TST"]['dnaSequence']))

pass