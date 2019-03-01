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

    def returnFirstCodingSequence(self, seq, maxPeptideLength):
        start = None
        peptideDNASeq = None
        for i in range(0, len(seq)):
            subSeq = seq[i:i+3]
            peptide = subSeq.translate(CodonTable.unambiguous_dna_by_id[1], to_stop=True)
            if str(peptide) == "M":
                start = i
                break
        startToEnd = seq[start:]
        endPeptide = None
        for i in range(0, len(startToEnd)//3):
            if maxPeptideLength is not None:
                if i+1 <= maxPeptideLength:
                    triplet = startToEnd[i*3:(i+1)*3]
                    translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1], to_stop=True)
                    if translatedTriplet == '*':
                        peptideDNASeq = startToEnd[:(i+1)*3]
                else:
                    peptideDNASeq = startToEnd[:i*3]
                    break
            else:
                triplet = startToEnd[i*3:(i+1)*3]
                translatedTriplet = triplet.translate(CodonTable.unambiguous_dna_by_id[1], to_stop=True)
                if translatedTriplet == '*':
                    peptideDNASeq = startToEnd[:(i+1)*3]

        return peptideDNASeq

    def deconstructImportedCDNASequence(self, sequence, sequenceIdentifier, maxPeptideLength=None):
        seq = Seq(sequence)
        codingSeq = self.returnFirstCodingSequence(seq, maxPeptideLength)
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

    def deconstructDNASequence(self, sequence, sequenceIdentifier, coding):
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

    def makeNewDeconstructedSequenceFromConstructSequenceWithPeptideMutation(self, deconstructedSeq, constructName, positionToMutate, beforeAminoAcid, afterAminoAcid):
        newDeconstructedList = []
        i = 0
        for item in deconstructedSeq['deconstructedList']:
            item['peptidePosition'] = i+1
            item['dnaStartPosition'] = i*3+1
            item['dnaEndPosition'] = (i+1)*3
            item['sequenceName'] = constructName
            if item['peptidePosition'] == positionToMutate:
                if item['peptide'] == beforeAminoAcid:
                    item['peptide'] = afterAminoAcid
                    item['dna'] = self.backTranslateDnaCodonMinimalNucelotideChanges(item['dna'], 'D')
            newDeconstructedList.append(item)
            i = i+1
        newDeconstructedSeq = {}
        newDeconstructedSeq['sequenceIdentifier'] = constructName
        seq = self.returnDnaSequenceFromDeconstructedList(newDeconstructedList)
        newDeconstructedSeq['dnaSequence'] = str(seq)
        newDeconstructedSeq['deconstructedList'] = newDeconstructedList
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.returnPeptideSequenceFromDeconstructedList(newDeconstructedSeq['deconstructedList'])
        self.allDeconstructedSequences[constructName] = newDeconstructedSeq
        self.allConstructs[constructName] = newDeconstructedSeq

    def makeNewDeconstructedSequenceFromDeconstructedSequencePeptideRange(self, deconstructedSeq, start, end, name):
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
        seq = self.returnDnaSequenceFromDeconstructedList(deconstructedRange)
        newDeconstructedSeq['dnaSequence'] = str(seq)
        newDeconstructedSeq['deconstructedList'] = deconstructedRange
        newDeconstructedSeq['coding'] = True
        newDeconstructedSeq['peptideSequence'] = self.returnPeptideSequenceFromDeconstructedList(newDeconstructedSeq['deconstructedList'])
        self.allDeconstructedSequences[name] = newDeconstructedSeq
        self.allConstructs[name] = newDeconstructedSeq

    def backTranslateDnaCodonMinimalNucelotideChanges(self, beforeDnaCodon, desiredAminoAcid):
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
                return codon





    def returnDnaSequenceFromDeconstructedList(self, deconstructedList):
        seq = Seq('')
        for item in deconstructedList:
            seq += item['dna']
        return seq

    def returnPeptideSequenceFromDeconstructedList(self, deconstructedList):
        seq = Seq('')
        for item in deconstructedList:
            if item['coding']:
                seq += item['peptide']
        return seq

    def returnPeptideSequenceRangeFromDeconstructedList(self, deconstructedList, start, end):
        seq = Seq('')
        for item in deconstructedList:
            if item['peptidePosition'] >= start and item['peptidePosition'] <= end:
                seq += Seq(item['peptide'])
        return seq

    def importSeqeuenceByNCBIIdentifier(self, ncbiSequenceIdentifier):
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

    def compareSequence(self, seq1, seq2):
        zipped = zip(seq1.lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1]
        return True

    def compareDNAConstructToSequence(self, construct, seq2):
        zipped = zip(construct['dnaSequence'].lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1], construct['deconstructedList'][i//3]['sequenceName']
        return True

    def comparePeptideConstructToSequence(self, construct, seq2):
        zipped = zip(construct['peptideSequence'].lower(), seq2.lower())
        zipList = list(zipped)
        for i in range(0, len(zipList)):
            if zipList[i][0] != zipList[i][1]:
                return False, i, zipList[i][0], zipList[i][1], construct['deconstructedList'][i]['sequenceName']
        return True

    def createSeqObjectFromString(self, seq):
        #remove whitespace
        seqWithoutWhitespace = ''.join(seq.split())
        seqObj = Seq(seqWithoutWhitespace)
        return seqObj.lower()

    def createConstructFromDeconstructedSequences(self, deconstructedSequenceNames, constructName):
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
        constructDict['dnaSequence']  = self.returnDnaSequenceFromDeconstructedList(constructDict['deconstructedList'])
        constructDict['peptideSequence'] = self.returnPeptideSequenceFromDeconstructedList(constructDict['deconstructedList'])
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