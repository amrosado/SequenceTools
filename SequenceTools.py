from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Data import CodonTable

class SequenceTools:
    currentSequence = None

    allSequences = None

    allDeconstructedSequences = None

    allEFetchRecords = None

    allConstructs = None

    def return_first_coding_sequence(self, seq, maxPeptideLength=None, minPeptideLength=None):
        start = None
        peptideDNASeq = None
        for i in range(0, len(seq)):
            subSeq = seq[i:i+3]
            peptide = subSeq.translate(CodonTable.unambiguous_dna_by_id[1], to_stop=True)
            seqPeptide = seq[i:].translate(CodonTable.unambiguous_dna_by_id[1], to_stop=True)
            if str(peptide) == "M":
                seqPeptideLen = len(seqPeptide)
                if minPeptideLength is not None and len(seqPeptide) >= minPeptideLength:
                    start = i
                    break
        startToEnd = seq[start:]
        endPeptide = None
        rang = range(len(startToEnd)//3)
        for i in range(len(startToEnd)//3):
            if maxPeptideLength is not None:
                if i+1 <= maxPeptideLength:
                    triplet = startToEnd[i*3:(i+1)*3]
                    translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1])
                    if translatedTriplet == '*':
                        peptideDNASeq = startToEnd[:(i+1)*3]
                else:
                    peptideDNASeq = startToEnd[:i*3]
                    break
            else:
                triplet = startToEnd[i*3:(i+1)*3]
                translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1])
                if translatedTriplet == '*':
                    peptideDNASeq = startToEnd[:(i+1)*3]
                    break
                elif i == (len(startToEnd)//3-1):
                    peptideDNASeq = startToEnd[:(i+1)*3]
        return peptideDNASeq

    def deconstruct_imported_cdna_sequence(self, sequence, sequenceIdentifier, maxPeptideLength=None, minPeptideLength=None):
        seq = Seq(sequence)
        codingSeq = self.return_first_coding_sequence(seq, maxPeptideLength, minPeptideLength)
        protein = codingSeq.translate(CodonTable.unambiguous_dna_by_id[1])
        deconstructedDict = {}
        deconstructedDict['coding'] = True
        deconstructedDict['sequenceIdentifier'] = sequenceIdentifier
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
            codonDict['sequenceName'] = sequenceIdentifier
            deconstructedDict['deconstructedList'].append(codonDict)
        self.allDeconstructedSequences[sequenceIdentifier] = deconstructedDict
        self.allConstructs[sequenceIdentifier] = deconstructedDict

    def deconstruct_dna_sequence(self, sequence, sequenceIdentifier, coding):
        if coding:
            protein = sequence.translate(CodonTable.unambiguous_dna_by_id[1])
            deconstructedDict = {}
            deconstructedDict['coding'] = coding
            deconstructedDict['sequenceIdentifier'] = sequenceIdentifier
            deconstructedDict['dnaSequence'] = sequence
            deconstructedDict['deconstructedList'] = []
            deconstructedDict['peptideSequence'] = Seq('')
            for i in range(0, len(sequence)//3):
                codonDict = {}
                codonDict['peptide'] = protein[i]
                codonDict['coding'] = coding
                deconstructedDict['peptideSequence'] += Seq(protein[i])
                codonDict['dna'] = sequence[i*3:(i+1)*3]
                codonDict['peptidePosition'] = i+1
                codonDict['dnaStartPosition'] = i*3+1
                codonDict['dnaEndPosition'] = (i+1)*3
                codonDict['sequenceName'] = sequenceIdentifier
                deconstructedDict['deconstructedList'].append(codonDict)
        else:
            deconstructedDict = {}
            deconstructedDict['coding'] = coding
            deconstructedDict['sequenceIdentifier'] = sequenceIdentifier
            deconstructedDict['dnaSequence'] = sequence
            deconstructedDict['deconstructedList'] = []
            for i in range(0, len(sequence)):
                codeDict = {}
                codeDict['dna'] = sequence[i]
                codeDict['coding'] = False
                codeDict['dnaStartPosition'] = i+1
                codeDict['dnaEndPosition'] = i+1
                codeDict['sequenceName'] = sequenceIdentifier
                deconstructedDict['deconstructedList'].append(codeDict)


        self.allDeconstructedSequences[sequenceIdentifier] = deconstructedDict

    def make_new_deconstructed_sequence_from_construct_sequence_with_peptide_mutation(self, deconstructedSeq, constructName, positionToMutate, beforeAminoAcid, afterAminoAcid):
        newDeconstructedList = []
        i = 0
        for item in deconstructedSeq['deconstructedList']:
            item['peptidePosition'] = i+1
            item['dnaStartPosition'] = i*3+1
            item['dnaEndPosition'] = (i+1)*3
            item['sequenceName'] = constructName
            if item['peptidePosition'] == positionToMutate:
                if item['peptide'] == beforeAminoAcid:
                    item['peptide'] = Seq(afterAminoAcid)
                    item['dna'] = self.back_translate_dna_codon_minimal_nucelotide_changes(item['dna'], 'D')
            newDeconstructedList.append(item)
            i = i+1
        newDeconstructedSeq = {}
        newDeconstructedSeq['sequenceIdentifier'] = constructName
        seq = self.return_dna_sequence_from_deconstructed_list(newDeconstructedList)
        newDeconstructedSeq['dnaSequence'] = str(seq)
        newDeconstructedSeq['deconstructedList'] = newDeconstructedList
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(newDeconstructedSeq['deconstructedList'])
        self.allDeconstructedSequences[constructName] = newDeconstructedSeq
        self.allConstructs[constructName] = newDeconstructedSeq

    def make_new_deconstructed_sequence_from_deconstructed_sequence_peptide_range(self, deconstructedSeq, start, end, name):
        deconstructedRange = []
        i = 0
        for item in deconstructedSeq['deconstructedList']:
            if item['peptidePosition'] >= start and item['peptidePosition'] <= end:
                item['peptidePosition'] = i+1
                item['dnaStartPosition'] = i*3+1
                item['dnaEndPosition'] = (i+1)*3
                item['sequenceName'] = name
                i+=1
                deconstructedRange.append(item)
        newDeconstructedSeq = {}
        newDeconstructedSeq['sequenceIdentifier'] = name
        seq = self.return_dna_sequence_from_deconstructed_list(deconstructedRange)
        newDeconstructedSeq['dnaSequence'] = str(seq)
        newDeconstructedSeq['deconstructedList'] = deconstructedRange
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(newDeconstructedSeq['deconstructedList'])
        self.allDeconstructedSequences[name] = newDeconstructedSeq
        self.allConstructs[name] = newDeconstructedSeq

    def back_translate_dna_codon_minimal_nucelotide_changes(self, beforeDnaCodon, desiredAminoAcid):
        standardTable = CodonTable.unambiguous_dna_by_id[1]
        possibleCodons = []
        codonScores = {}
        bestCodonScore = 0
        for forwardCodon in standardTable.forward_table:
            if standardTable.forward_table[forwardCodon] == desiredAminoAcid:
                possibleCodons.append(forwardCodon)
        for possibleCodon in possibleCodons:
            zipped = zip(beforeDnaCodon.upper(), possibleCodon)
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
        return seq

    def return_peptide_sequence_from_deconstructed_list(self, deconstructedList):
        seq = Seq('')
        for item in deconstructedList:
            if item['coding']:
                seq += item['peptide']
        return seq

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
                eFetchHandle = Entrez.efetch(db="nucleotide", id=id, retmode="xml")
                eFetchRecord = Entrez.read(eFetchHandle)
                eFetchHandle.close()
                for result in eFetchRecord:
                    self.allEFetchRecords[ncbiSequenceIdentifier] = result
                    self.allSequences[ncbiSequenceIdentifier] = result['GBSeq_sequence']

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
            deconstuctedSeq = self.allDeconstructedSequences[name]
            constructDict['deconstructedComponents'].append(deconstuctedSeq)
            if self.allDeconstructedSequences[name]['coding']:
                for item in deconstuctedSeq['deconstructedList']:
                    item['peptidePosition'] = i/3+1
                    item['dnaStartPosition'] = i+1
                    item['dnaEndPosition'] = i+3
                    i+=3
                    constructDict['deconstructedList'].append(item)
            else:
                for item in deconstuctedSeq['deconstructedList']:
                    item['dnaStartPosition'] = i
                    item['dnaEndPosition'] = i
                    i+=1
                    constructDict['deconstructedList'].append(item)
        constructDict['dnaSequence']  = self.return_dna_sequence_from_deconstructed_list(constructDict['deconstructedList'])
        constructDict['peptideSequence'] = self.return_peptide_sequence_from_deconstructed_list(constructDict['deconstructedList'])
        self.allConstructs[constructName] = constructDict


    def __init__(self, apiKey=None, email=None):
        Entrez.tool = "SequenceTools"
        if apiKey:
            Entrez.api_key = apiKey
        if email:
            Entrez.email = email

        self.allSequences = {}
        self.allESearchRecords = {}
        self.allEFetchRecords = {}
        self.allDeconstructedSequences = {}
        self.allConstructs = {}