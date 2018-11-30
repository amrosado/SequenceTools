# Protein Construct Tools

## Description
A small python library written to import cDNA sequences from NCBI database and form chimeric protein constructs.

## Dependencies
* Biopython

## Install
* Written for python 3

* Install requirements
    > pip install -r requirements.txt

## Use
1) Import Sequence Tools

    `from SequenceTools import SequenceTools`

2) Create a SequenceTools object with your email

    `tools = SequenceTools(email="your@email.com")`
    
    For NCBI api key:
    
    `tools = SequenceTools(apiKey="api_key", email="your@email.com")`

3) Import Sequence from NCBI UniGene

    `tools.importSeqeuenceByNCBIIdentifier("AK079513.1")`

4) Deconstruct imported sequence and give deconstruction key

    `tools.deconstructImportedCDNASequence(tools.allSequences["AK079513.1"], "AK079513.1")`

5) Make new deconstructed sequence from range to isolate portion (in this example the extracellular peptides) and give key.

    `tools.makeNewDeconstructedSequenceFromDeconstructedSequencePeptideRange(tools.allDeconstructedSequences["AK079513.1"], 24, 244, "CD86_Extracellular")`

6) Create sequences from copied DNA of other proteins, etc.

    `nhe1Seq = tools.createSeqObjectFromString("GCTAGC")`

7) Deconstruct sequences, give key for reference, and specify if coding

    `tools.deconstructDNASequence(nhe1Seq, "NheI", False)`
 
8) Make a new construct from deconstructed sequences and give a key

    `tools.createConstructFromDeconstructedSequences(['NheI' ,'SecretionSignal', 'Linker1', 'CD86_Extracellular', "Linker2", 'APTag', 'Linker3', 'TEV', 'HIS','STOPS', 'EcoRI'], 'NheI-CD86Extracellular-APTag-TEV-HIS-EcoRI')`
    
9) Compare construct with seq

    `comp5 = tools.compareDNAConstructToSequence(tools.allConstructs['NheI-CD86Extracellular-APTag-TEV-HIS-EcoRI'], genScriptSeq)`    

## Contact
For new feature requests or debugging contact: arosad2@protonmail.ch.

