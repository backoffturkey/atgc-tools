# tcl blast parser usage

# TCL Blast Parser Usage #

Note about NCBI Blast usage: it's important to run BLAST search with options:

**"-V T"**: -V  Force use of the legacy BLAST engine (T/F); Optional; default = F

**"-C F"**: Use (false) composition-based statistics for blastp or tblastn

for example:

  * Blast-X search:
blastall -p blastx -V T -F F -e 1e-20 -b 24 -v 24 -d db\_file -i input.fasta -o BlastX.out

  * Blast-P search:
blastall -p blastp -C F -V T -F F -e 1e-20 -b 24 -v 24 -d db\_file -i input.fasta -o BlastP.out

  * Blast-N search:
blastall -p blastn -V T -F F -e 1e-20 -b 24 -v 24 -d db\_file -i input.fasta -o BlastN.out _(for long DNA sequences)_

blastall -p blastn -V T -F F -e 0.001 -b 10000 -v 10000 -d db\_file -i input.fasta -o BlastN.out _(for short DNA fragments)_

**tcl BLAST parser run options:**

tclsh tcl\_blast\_parser\_123\_V047.tcl BlastX.out BlastX.out.parsed 20 40 100 MATRIX

  * BlastX.out - blast output file
  * BlastX.out.parsed - file name prefix for parsed output files
  * 20 - expect cutoff for MATRIX files
  * 40 - identity cutoff for MATRIX files
  * 100 - alignment length cutoff for MATRIX files

# Tcl Blast Parser Output Files #

**"`*`.all\_hits" file (tab-delimited table) contains info about ALL HITS in blast report:**

  * 1 column: "Query" sequence ID
  * 2 column: "Subject" sequence ID
  * 3 column: normalized expectation (-log(Exp))
  * 4 column: percent of identity
  * 5 column: number of perfect matches
  * 6 column: length of the alignment
  * 7 column: hit number for primary alignment (1 is the best first hit)
  * 8 column: hit number for alternative alignment
  * 9 column: PRM stands for primary alignment, ALT for alternative alignment
  * 10 column: frame of translation (+1, +2, +3, -1, -2 or -3)
  * 11 column: first position of the "Query" sequence in the alignment
  * 12 column: last position of the "Query" sequence in the alignment
  * 13 column: first position of the "Subject" sequence in the alignment
  * 14 column: last position of the "Subject" sequence in the alignment
  * 15 column: "`***`" - visual mark
  * 16 column: "Query" length
  * 17 column: "Subject" length
  * 18 column: Gap info (GAPPED or NO\_GAPS)


**"`*`.info1" file (tab-delimited table) contains info about FIRST BEST hits only and description line of the subject from BLAST report:**

  * 1 column: "Query" sequence ID
  * 2 column: "Subject" sequence ID
  * 3 column: description line for "Subject" sequence
  * 4 column: normalized expectation (-log(Exp)/100) (note, that expectation is normalized between 0 and 1)
  * 5 column: percent of identity
  * 6 column: number of perfect matches
  * 7 column: length of the alignment
  * 8 column: hit number for primary alignment (1 is the best first hit)
  * 9 column: "Query" length
  * 10 column: "Subject" length


**"`*`.info2" file (tab-delimited table) contains info about ALL PRIMARY HITS in blast report. It is almost identical to All Hits file with exception that it contains description line for subject, expectation is normalized between 0 and 1 and there is no information for ALT - alternative alignments:**

  * 1 column: "Query" sequence ID
  * 2 column: "Subject" sequence ID
  * 3 column: description line for "Subject" sequence
  * 4 column: normalized expectation (-log(Exp)/100) (note, that expectation is normalized between 0 and 1)
  * 5 column: percent of identity
  * 6 column: number of perfect matches
  * 7 column: length of the alignment
  * 8 column: hit number for primary alignment (1 is the best first hit)
  * 9 column: frame of translation (+1, +2, +3, -1, -2 or -3)
  * 10 column: first position of the "Query" sequence in the alignment
  * 11 column: last position of the "Query" sequence in the alignment
  * 12 column: first position of the "Subject" sequence in the alignment
  * 13 column: last position of the "Subject" sequence in the alignment
  * 14 column: Length of "Query"(nt or aa) / Length of "Subject"(nt or aa)
  * 15 column: Gap info (GAPPED or NO\_GAPS)


**MATRIX output files (expect and identity):**

**"`*`.matrix.expect":**

  * 1st and 2nd column: pair of sequence IDs
  * 3rd column: EXPECT value normalized between 0 and 1 (-log(Exp))/100
  * 4th column: identity values (0-100 range)
  * 5th column: length of the alignment

**"`*`.matrix.identity":**

  * 1st and 2nd column: pair of sequence IDs
  * 3rd column: IDENTITY value normalized between 0 and 1 (identity`%` / 100)
  * 4th column: normalized expectation (-log(Exp)) (0-100 range)
  * 5th column: length of the alignment