from SequenceTools import SequenceTools

class CommonSeqsTags:
    def __init__(self):
        self.seq_tools = SequenceTools()

    def return_spe1(self):
        spe_1 = self.seq_tools.create_seq_object_from_string("ACTAGT")

        self.seq_tools.deconstruct_dna_sequence(spe_1, "SpeI", False)

        return self.seq_tools.return_deconstructed_sequence("SpeI")

    def return_xba1(self):
        xba_1 = self.seq_tools.create_seq_object_from_string("TCTAGA")

        self.seq_tools.deconstruct_dna_sequence(xba_1, "XbaI", False)

        return self.seq_tools.return_deconstructed_sequence("XbaI")

    def return_xho1(self):
        xho_1 = self.seq_tools.create_seq_object_from_string("CTCGAG")

        self.seq_tools.deconstruct_dna_sequence(xho_1, "XhoI", False)

        return self.seq_tools.return_deconstructed_sequence("XhoI")

    def return_kozak_non_coding(self):
        kozak_seq_non_coding = self.seq_tools.create_seq_object_from_string("gccAcc")

        self.seq_tools.deconstruct_dna_sequence(kozak_seq_non_coding, "KozakNonCoding", False)

        return self.seq_tools.return_deconstructed_sequence("KozakNonCoding")

    def return_kozak_coding(self):
        kozak_seq_coding = self.seq_tools.create_seq_object_from_string("ATGGGT")

        self.seq_tools.deconstruct_dna_sequence(kozak_seq_coding, "KozakCoding", True)

        return self.seq_tools.return_deconstructed_sequence("KozakCoding")

    def return_ssmigk(self):
        # Secretion signal mouse Ig kappa chain
        ssmigk = self.seq_tools.create_seq_object_from_string("ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGAC")

        self.seq_tools.deconstruct_dna_sequence(ssmigk, "SSMIGGK", True)

        return self.seq_tools.return_deconstructed_sequence("SSMIGGK")

    def return_ap_tag(self):
        ap = self.seq_tools.create_seq_object_from_string("GGTCTGAATGATATTTTCGAAGCGCAGAAAATTGAATGGCATGAA")

        self.seq_tools.deconstruct_dna_sequence(ap, "APTag", True)

        return self.seq_tools.return_deconstructed_sequence("APTag")

    def return_tev(self):
        tev = self.seq_tools.create_seq_object_from_string("GAGAACCTATACTTCCAAGGA")

        tev_seq = self.seq_tools.deconstruct_imported_cdna_sequence(tev, "TEV", True)

        return self.seq_tools.return_deconstructed_sequence("TEV")

    def return_tst_tag(self):
        tst = self.seq_tools.create_seq_object_from_string("TGGTCCCACCCCCAGTTCGAGAAGGGCGGCGGTAGTGGAGGGGGCGGATCTGGCGGCTCAGCTTGGAGCCACCCCCAGTTCGAAAAG")

        self.seq_tools.deconstruct_dna_sequence(tst, "TST", True)

        return self.seq_tools.return_deconstructed_sequence("TST")

    def return_8xhis_tag(self):
        his8x = self.seq_tools.create_seq_object_from_string("CATCACCACCATCACCACCATCAT")

        self.seq_tools.deconstruct_dna_sequence(his8x, "8xHis", True)

        return self.seq_tools.return_deconstructed_sequence("8xHis")

    def return_stops(self):
        stops = self.seq_tools.create_seq_object_from_string("TAGTAA")

        self.seq_tools.deconstruct_dna_sequence(stops, "Stops", True)

        return self.seq_tools.return_deconstructed_sequence("Stops")

    def return_higg1_fc(self):
        peptide_seq = """
        ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSS
        GLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGG
        PSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYN
        STYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDE
        LTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRW
        QQGNVFSCSVMHEALHNHYTQKSLSLSPELQLEESCAEAQDGELDGLWTTITIFITLFLL
        SVCYSATVTFFKVKWIFSSVVDLKQTIIPDYRNMIGQGA
        """

        peptide_seq = peptide_seq.replace("\n", "").replace(" ", "")

        cdna_seq = self.seq_tools.back_translate_amino_acid_sequence(peptide_seq)

        self.seq_tools.deconstruct_dna_sequence(cdna_seq, "hIgG1", True)

        self.seq_tools.make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range("hIgG1", 121, 325, "hIgG1_Fc")

        return self.seq_tools.return_deconstructed_sequence("hIgG1_Fc")

