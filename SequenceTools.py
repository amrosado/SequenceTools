from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Sequencing import Ace

import copy
import random
import io

class SequenceTools:
    current_sequence = None

    all_sequences = None

    all_deconstructed_sequences = None

    all_efetch_records = None

    all_constructs = None

    def import_deconstructed_sequence(self, deconstructed_seq):
        self.all_deconstructed_sequences[deconstructed_seq['sequenceIdentifier']] = deconstructed_seq

    def return_first_coding_sequence(self, seq, max_peptide_length=None, min_peptide_length=None):
        start = None
        peptide_dna_seq = None

        for i in range(0, len(seq)):
            sub_seq = seq[i:i+3]
            peptide = sub_seq.translate(table='Standard', to_stop=True)
            seq_peptide = seq[i:].translate(table='Standard', to_stop=True)
            if str(peptide) == "M":
                if min_peptide_length is not None and len(seq_peptide) >= min_peptide_length:
                    start = i
                    break
        start_to_end = seq[start:]
        end_peptide = None
        rang = range(len(start_to_end)//3)
        for i in range(len(start_to_end)//3):
            if max_peptide_length is not None:
                if i+1 <= max_peptide_length:
                    triplet = start_to_end[i*3:(i+1)*3]
                    translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1])
                    if translatedTriplet == '*':
                        peptide_dna_seq = start_to_end[:(i+1)*3]
                else:
                    peptide_dna_seq = start_to_end[:i*3]
                    break
            else:
                triplet = start_to_end[i*3:(i+1)*3]
                translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1])
                if translatedTriplet == '*':
                    peptide_dna_seq = start_to_end[:(i+1)*3]
                    break
                elif i == (len(start_to_end)//3-1):
                    peptide_dna_seq = start_to_end[:(i+1)*3]
        return peptide_dna_seq

    def deconstruct_imported_orf_sequence(self, sequence, sequence_identifier, sequence_to_search_for, min_peptide_length=None):
        seq = Seq(sequence)
        for i in range(0, len(seq)):
            potential_protein = seq[i:].translate(CodonTable.unambiguous_dna_by_id[1])
            find_result = str(potential_protein.lower()).find(sequence_to_search_for.lower())
            if find_result > 0:
                coding_seq = self.return_first_coding_sequence(seq[find_result*3:], min_peptide_length=min_peptide_length)
                protein = coding_seq.translate(CodonTable.unambiguous_dna_by_id[1])
                deconstructedDict = {}
                deconstructedDict['coding'] = True
                deconstructedDict['sequenceIdentifier'] = sequence_identifier
                deconstructedDict['dnaSequence'] = coding_seq
                deconstructedDict['peptideSequence'] = Seq('')
                deconstructedDict['deconstructedList'] = []
                for i in range(0, len(coding_seq)//3):
                    codonDict = {}
                    codonDict['peptide'] = protein[i]
                    deconstructedDict['peptideSequence'] += Seq(protein[i])
                    codonDict['dna'] = coding_seq[i*3:(i+1)*3]
                    codonDict['peptidePosition'] = i+1
                    codonDict['dnaStartPosition'] = i*3+1
                    codonDict['dnaEndPosition'] = (i+1)*3
                    codonDict['coding'] = True
                    codonDict['sequenceName'] = sequence_identifier
                    deconstructedDict['deconstructedList'].append(codonDict)
                self.all_deconstructed_sequences[sequence_identifier] = deconstructedDict
                self.all_constructs[sequence_identifier] = deconstructedDict
                break


    def deconstruct_imported_cdna_sequence(self, sequence, sequence_identifier, max_peptide_length=None, min_peptide_length=None):
        seq = Seq(sequence)
        codingSeq = self.return_first_coding_sequence(seq, max_peptide_length, min_peptide_length)
        protein = codingSeq.translate(CodonTable.unambiguous_dna_by_id[1])
        deconstructedDict = {}
        deconstructedDict['coding'] = True
        deconstructedDict['sequenceIdentifier'] = sequence_identifier
        deconstructedDict['dnaSequence'] = codingSeq
        deconstructedDict['peptideSequence'] = Seq('')
        deconstructedDict['deconstructedList'] = []
        for i in range(0, len(codingSeq)//3):
            codonDict = {}
            codonDict['peptide'] = protein[i]
            deconstructedDict['peptideSequence'] += Seq(protein[i])
            codonDict['dna'] = codingSeq[i*3:(i+1)*3]
            codonDict['peptidePosition'] = i+1
            codonDict['dnaStartPosition'] = i*3+1
            codonDict['dnaEndPosition'] = (i+1)*3
            codonDict['coding'] = True
            codonDict['sequenceName'] = sequence_identifier
            deconstructedDict['deconstructedList'].append(codonDict)
        self.all_deconstructed_sequences[sequence_identifier] = deconstructedDict
        self.all_constructs[sequence_identifier] = deconstructedDict

    def deconstruct_dna_sequence(self, sequence, sequence_identifier, coding):
        if coding:
            protein = sequence.translate(CodonTable.unambiguous_dna_by_id[1])
            deconstructedDict = {}
            deconstructedDict['coding'] = coding
            deconstructedDict['sequenceIdentifier'] = sequence_identifier
            deconstructedDict['dnaSequence'] = sequence
            deconstructedDict['deconstructedList'] = []
            deconstructedDict['peptideSequence'] = Seq('')
            for i in range(0, len(sequence)//3):
                codonDict = {}
                codonDict['peptide'] = protein[i].lower()
                codonDict['coding'] = coding
                deconstructedDict['peptideSequence'] += Seq(protein[i])
                codonDict['dna'] = sequence[i*3:(i+1)*3].lower()
                codonDict['peptidePosition'] = i+1
                codonDict['dnaStartPosition'] = i*3+1
                codonDict['dnaEndPosition'] = (i+1)*3
                codonDict['sequenceName'] = sequence_identifier
                deconstructedDict['deconstructedList'].append(codonDict)
        else:
            deconstructedDict = {}
            deconstructedDict['coding'] = coding
            deconstructedDict['sequenceIdentifier'] = sequence_identifier
            deconstructedDict['dnaSequence'] = sequence
            deconstructedDict['deconstructedList'] = []
            for i in range(0, len(sequence)):
                codeDict = {}
                codeDict['dna'] = sequence[i].lower()
                codeDict['coding'] = False
                codeDict['dnaStartPosition'] = i+1
                codeDict['dnaEndPosition'] = i+1
                codeDict['sequenceName'] = sequence_identifier
                deconstructedDict['deconstructedList'].append(codeDict)

        self.all_constructs[sequence_identifier] = deconstructedDict
        self.all_deconstructed_sequences[sequence_identifier] = deconstructedDict

    def make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(self, deconstructedSeq, constructName, positionToMutate, beforeAminoAcid, afterAminoAcid):
        newDeconstructedList = []
        i = 0

        for item in deconstructedSeq['deconstructedList']:
            item_copy = item.copy()

            item_copy['peptidePosition'] = i+1
            item_copy['dnaStartPosition'] = i*3+1
            item_copy['dnaEndPosition'] = (i+1)*3
            item_copy['sequenceName'] = constructName
            if item_copy['peptidePosition'] == positionToMutate:
                if item['peptide'] == beforeAminoAcid:
                    item_copy['peptide'] = Seq(afterAminoAcid)
                    item_copy['dna'] = self.back_translate_dna_codon_minimal_nucelotide_changes(item_copy['dna'], afterAminoAcid).lower()

            newDeconstructedList.append(item_copy)
            i = i+1

        newDeconstructedSeq = {}
        newDeconstructedSeq['sequenceIdentifier'] = constructName
        seq = self.return_dna_sequence_from_deconstructed_list(newDeconstructedList)
        newDeconstructedSeq['dnaSequence'] = str(seq)
        newDeconstructedSeq['deconstructedList'] = newDeconstructedList
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(newDeconstructedSeq['deconstructedList'])

        self.all_deconstructed_sequences[constructName] = newDeconstructedSeq
        self.all_constructs[constructName] = newDeconstructedSeq

    def make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(self, deconstructed_seq, start, end, name):
        if isinstance(deconstructed_seq, str):
            deconstructed_seq = self.all_deconstructed_sequences[deconstructed_seq]

        deconstructedRange = []
        i = 0

        for item in deconstructed_seq['deconstructedList']:
            if item['peptidePosition'] >= start and item['peptidePosition'] <= end:
                item_copy = item.copy()

                item_copy['peptidePosition'] = i+1 - start
                item_copy['dnaStartPosition'] = (i+1)*3 - start*3
                item_copy['dnaEndPosition'] = (i+2)*3 - start*3
                item_copy['sequenceName'] = name

                deconstructedRange.append(item_copy)

            i+=1



        newDeconstructedSeq = {}
        newDeconstructedSeq['sequenceIdentifier'] = name
        seq = self.return_dna_sequence_from_deconstructed_list(deconstructedRange)
        newDeconstructedSeq['dnaSequence'] = str(seq).lower()
        newDeconstructedSeq['deconstructedList'] = deconstructedRange
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(newDeconstructedSeq['deconstructedList'])

        self.all_deconstructed_sequences[name] = newDeconstructedSeq
        self.all_constructs[name] = newDeconstructedSeq

    def back_translate_amino_acid_sequence(self, desired_amino_acid_seq):
        dna_seq = Seq('')

        for i in desired_amino_acid_seq:
            dna_seq += self.back_translate_amino_acid_to_codon(i)

        return dna_seq

    def back_translate_amino_acid_sequence_random_codon(self, desired_amino_acid_seq):
        dna_seq = Seq('')

        for i in desired_amino_acid_seq:
            dna_seq += self.back_translate_amino_acid_to_random_codon(i)

        return dna_seq

    def back_translate_amino_acid_to_random_codon(self, desired_amino_acid):
        standard_table = CodonTable.unambiguous_dna_by_id[1]

        potential_codon_list = []

        for codon in standard_table.forward_table:
            if standard_table.forward_table[codon] == desired_amino_acid.upper():
                potential_codon_list.append(codon)

        random_int = random.randint(0, len(potential_codon_list) - 1)

        return potential_codon_list[random_int]

    def back_translate_amino_acid_to_codon(self, desired_amino_acid):
        standard_table = CodonTable.unambiguous_dna_by_id[1]

        return standard_table.back_table[desired_amino_acid.upper()]

    def back_translate_dna_codon_minimal_nucelotide_changes(self, before_dna_codon, desired_amino_acid):
        standardTable = CodonTable.unambiguous_dna_by_id[1]
        possibleCodons = []
        codonScores = {}
        bestCodonScore = 0

        for forwardCodon in standardTable.forward_table:
            if standardTable.forward_table[forwardCodon] == desired_amino_acid:
                possibleCodons.append(forwardCodon)

        for possibleCodon in possibleCodons:
            zipped = zip(before_dna_codon.upper(), possibleCodon)
            score = 0
            for zippedItem in zipped:
                if zippedItem[0] == zippedItem[1]:
                    score+=1
            codonScores[possibleCodon] = score

        bestCodonScore = max(codonScores.values())

        for codon in codonScores:
            if bestCodonScore == codonScores[codon]:
                return Seq(codon)

    def return_dna_sequence_from_deconstructed_list(self, deconstructedList):
        seq = Seq('')
        for item in deconstructedList:
            seq += item['dna']
        return seq.lower()

    def return_peptide_sequence_from_deconstructed_list(self, deconstructedList):
        seq = Seq('')
        for item in deconstructedList:
            if item['coding']:
                seq += item['peptide']
        return seq.lower()

    def insert_sequence_into_peptide_position(self, deconstructed_list_name, insert_deconstructed_list_name, new_list_name, insert_start_position):
        list_to_insert_into = copy.deepcopy(self.all_deconstructed_sequences[deconstructed_list_name])
        insert_list_reverse = copy.deepcopy(self.all_deconstructed_sequences[insert_deconstructed_list_name])
        insert_list_reverse['deconstructedList'].reverse()

        for j in insert_list_reverse['deconstructedList']:
            list_to_insert_into['deconstructedList'].insert(insert_start_position, j)

        dna_seq = self.return_dna_sequence_from_deconstructed_list(list_to_insert_into['deconstructedList'])

        self.deconstruct_dna_sequence(dna_seq, new_list_name, True)

    def return_peptide_sequence_range_from_deconstructed_list(self, deconstructedList, start, end):
        seq = Seq('')
        for item in deconstructedList:
            if item['peptidePosition'] >= start and item['peptidePosition'] <= end:
                seq += Seq(item['peptide'])
        return seq

    def import_sequence_by_ncbi_identifier(self, ncbiSequenceIdentifier):
        if isinstance(ncbiSequenceIdentifier, str):
            eSearchHandle = Entrez.esearch(db="nucleotide", term=ncbiSequenceIdentifier)
            eSearchRecord = Entrez.read(eSearchHandle)
            eSearchHandle.close()
            for id in eSearchRecord['IdList']:
                eFetchHandle = Entrez.efetch(db="nucleotide", id=id, retype="gb", retmode="xml", strand="true")
                eFetchRecord = Entrez.read(eFetchHandle, validate=False)

                for result in eFetchRecord:
                    self.all_efetch_records[ncbiSequenceIdentifier] = result
                    if 'GBSeq_sequence' in result:
                        self.all_sequences[ncbiSequenceIdentifier] = result['GBSeq_sequence']
                    else:
                        with io.StringIO() as f:
                            f.write(result['GBSeq_contig'])
                            f.seek(0)
                            print(f.read())
                            gen = Ace.parse(f)
                            contig = next(gen)
                            pass

                eFetchHandle.close()

    def compare_sequence(self, seq1, seq2):
        zipped = zip(seq1.lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1]
        return True

    def compare_dna_construct_to_sequence(self, construct, seq2):
        zipped = zip(construct['dnaSequence'].lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1], construct['deconstructedList'][i//3]['sequenceName']
        return True

    def compare_peptide_construct_to_sequence(self, construct, seq2):
        if isinstance(construct, str):
            construct = self.all_deconstructed_sequences[construct]

        zipped = zip(construct['peptideSequence'].lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1], construct['deconstructedList'][i]['sequenceName']
        return True

    def create_seq_object_from_string(self, seq):
        #remove whitespace
        seqWithoutWhitespace = ''.join(seq.split())
        seqObj = Seq(seqWithoutWhitespace)
        return seqObj.lower()

    def create_construct_from_deconstructed_sequences(self, deconstructedSequenceNames, constructName):
        i = 0
        constructDict = {}
        constructDict['name'] = constructName
        constructDict['deconstructedComponentNames'] = deconstructedSequenceNames
        constructDict['deconstructedComponents'] = []
        constructDict['deconstructedList'] = []
        for name in deconstructedSequenceNames:
            deconstuctedSeq = self.all_deconstructed_sequences[name]
            constructDict['deconstructedComponents'].append(deconstuctedSeq)
            if self.all_deconstructed_sequences[name]['coding']:
                for item in deconstuctedSeq['deconstructedList']:
                    item_copy = item.copy()
                    item_copy['peptidePosition'] = i/3+1
                    item_copy['dnaStartPosition'] = i+1
                    item_copy['dnaEndPosition'] = i+3
                    i+=3
                    constructDict['deconstructedList'].append(item_copy)
            else:
                for item in deconstuctedSeq['deconstructedList']:
                    item_copy = item.copy()
                    item_copy['dnaStartPosition'] = i
                    item_copy['dnaEndPosition'] = i
                    i+=1
                    constructDict['deconstructedList'].append(item_copy)

        constructDict['dnaSequence']  = self.return_dna_sequence_from_deconstructed_list(constructDict['deconstructedList'])
        constructDict['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(constructDict['deconstructedList'])
        self.all_constructs[constructName] = constructDict

    def return_deconstructed_sequence(self, name):
        return self.all_deconstructed_sequences[name]


    def __init__(self, apiKey=None, email=None):
        Entrez.tool = "SequenceTools"
        if apiKey:
            Entrez.api_key = apiKey
        if email:
            Entrez.email = email

        self.all_sequences = {}
        self.allESearchRecords = {}
        self.all_efetch_records = {}
        self.all_deconstructed_sequences = {}
        self.all_constructs = {}