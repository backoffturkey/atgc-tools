#seqs\_processor\_and\_translator usage

# seqs\_processor\_and\_translator usage #

**seqs\_processor\_and\_translator\_bin\_V136\_AGCT.py version:**

```
bash-2.03$ python seqs_processor_and_translator_bin_V136_AGCT.py

  Program usage:
  input_file  output_file  DNA/prot  trans_frame[0 1 3 6]  genetic_code[1 2 3]  BIN/NOBIN  seqs_min_len SEQS_SPLIT/SEQS
  Script counts "ATGC" content in FASTA file
  and translate DNA sequence into protein
  0 - no translation; 1 - first frame
  3 - three frames; 6 - all six frame
  Genetic Code: 1 - Standard; 2 - Vertebrate Mitochondrial; 3 - Yeast Mitochondrial
  BIN option - to generate 0100111010101001100 style DNA file
  SEQS_SPLIT option - to split FASTA file into set of individual files/sequences

  use argument "help" - for help

bash-2.03$
```

```
   SEQS PROCESSOR HELP:
   -------------------

   Input file must be in FASTA format. Description of output:

   if no translation frame is chosen "0" then four output files will be generated:
   1. *.fasta - contains processed sequences in FASTA format
   2. *.stat  - contains information about GC% content per sequence
   3. *.tab   - tab delimited file with sequences
   4. *.01BIN - binary DNA representation (A-00 T-11 G-01 C-10)

   if first translation frame is chosen "1" then eleven output files will be generated:
   1. *.fasta - contains processed sequences in FASTA format
   2. *.forwrd - forward orientation, all non ATGCN letters replaced by "N"
   3. *.revcom - reverse-compliment orientation based on "forward" sequences
   4. *.stat  - contains information about GC% content per sequence
   5. *.tab   - tab delimited file with sequences
   6. *.tr_frame1 - translated sequences - frame 1
   7. *.tr_longest_frame - extracted longest frame
   8. *.tr_single_frame - extracted sequences with single ORF
   9. *.x0_log - log file with general info and error messages
  10. *.x1_log - log file with statistical info about translation - frame 1
  11. *.01BIN - binary DNA representation (A-00 T-11 G-01 C-10)

   if three translation frames are chosen "3" then sixteen output files will be generated:
   1. *.fasta - contains processed sequences in FASTA format
   2. *.forwrd - forward orientation, all non ATGCN letters replaced by "N"
   3. *.revcom - reverse-compliment orientation based on "forward" sequences
   4. *.stat  - contains information about GC% content per sequence
   5. *.tab   - tab delimited file with sequences
   6. *.tr_frame1 - translated sequences - frame 1
   7. *.tr_frame2 - translated sequences - frame 2
   8. *.tr_frame3 - translated sequences - frame 3
   9. *.tr_longest_frame - extracted longest frame
  10. *.tr_single_frame - extracted sequences with single ORF
  11. *.x0_log - log file with general info and error messages
  12. *.x1_log - log file with statistical info about translation - frame 1
  13. *.x2_log - log file with statistical info about translation - frame 2
  14. *.x3_log - log file with statistical info about translation - frame 3
  15. *.xx_log - log file with statistical info about all three translation ORFs
  16. *.01BIN - binary DNA representation (A-00 T-11 G-01 C-10)

   if six translation frames are chosen "6" then twenty two output files will be generated:
   1. *.fasta - contains processed sequences in FASTA format
   2. *.forwrd - forward orientation, all non ATGCN letters replaced by "N"
   3. *.revcom - reverse-compliment orientation based on "forward" sequences
   4. *.stat  - contains information about GC% content per sequence
   5. *.tab   - tab delimited file with sequences
   6. *.tr_frame1 - translated sequences - frame 1
   7. *.tr_frame2 - translated sequences - frame 2
   8. *.tr_frame3 - translated sequences - frame 3
   9. *.tr_frame4 - translated sequences - frame 4
  10. *.tr_frame5 - translated sequences - frame 5
  11. *.tr_frame6 - translated sequences - frame 6
  12. *.tr_longest_frame - extracted longest frame
  13. *.tr_single_frame - extracted sequences with single ORF
  14. *.x0_log - log file with general info and error messages
  15. *.x1_log - log file with statistical info about translation - frame 1
  16. *.x2_log - log file with statistical info about translation - frame 2
  17. *.x3_log - log file with statistical info about translation - frame 3
  18. *.x4_log - log file with statistical info about translation - frame 4
  19. *.x5_log - log file with statistical info about translation - frame 5
  20. *.x6_log - log file with statistical info about translation - frame 6
  21. *.xx_log - log file with statistical info about all six translation ORFs
  22. *.01BIN - binary DNA representation (A-00 T-11 G-01 C-10)

  Note, log files can be viewed using Excel like spreadsheet editor
```

**sequence processor changes for seqs\_processor\_and\_translator\_bin\_V144\_AGCT.py version:**

http://code.google.com/p/atgc-tools/source/browse/#svn/trunk

1. fixed bug for processing of short sequences (script went into infinite loop if last sequence shorter than cutoff value)

2. optional check for duplicated IDs (test for redundant IDs), for large sets (more than 100,000 sequences) it's reasonable to run without redundancy check, it works faster.

3. there is a new option/argument: RDN\_Y or RDN\_N - ID redundancy check (see above)

example run:

```
python seqs_processor_and_translator_bin_V144_AGCT.py file_in file_out DNA 0 1 NOBIN 24 SEQS RDN_Y
```
(slow with check for redundant IDs)

```
python seqs_processor_and_translator_bin_V144_AGCT.py file_in file_out DNA 0 1 NOBIN 24 SEQS RDN_N
```
(fast without test for duplicates in IDs)
