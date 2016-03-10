# Sequence Subgroup Extractor #

sequence\_subgroup\_extractor **seqs\_subgroup\_extr\_012.py** script works with two alternative options (FULL or RANGE) to extract subgroup of sequences from FASTA file:
  * sequences listed in a group file according to their IDs (option **FULL**)
  * fragments of sequences listed in a group file according to their IDs and length range (option **RANGE**)

by running the script without options it displays order of parameters/options to use:

```
bash-2.03$ python seqs_subgroup_extr_012.py

Program usage:
[input_file] [output_file] [group_id_file] [seqs_range]

Sequences listed in group_file will be extracted according to FASTA header IDs

bash-2.03$
```

**input\_file** - input file with sequences in FASTA format

**output\_file** - output file name

**group\_id\_file** - file with list of IDs to extract (and optional length range)

**seqs\_range** - two possible options: **FULL** or **RANGE**

Option **FULL** - sequence IDs from first column in a group file will be used to extract matched IDs and corresponding full-length sequences from input FASTA file, for example:

```
bash-2.03$ python seqs_subgroup_extr_012.py input_file output_file group_id_file FULL
```

Option **RANGE** - sequence IDs from first column in a group file will be used to extract matched IDs and corresponding sequence fragments from input FASTA file (sequence length range is given in second and third columns of group\_id\_file), for example:

```
bash-2.03$ python seqs_subgroup_extr_012.py input_file output_file group_id_file RANGE
```

group ID file for option FULL can be just a list of IDs:

```
ID_1
ID_2
ID_3
```

group ID file for option RANGE should be tab-delimited file with coordinates of sequence fragments to extract, for example:

```
ID_1    0   112
ID_2   11   200
ID_3   24   300
```