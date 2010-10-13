#!/usr/bin/python

####################################################
#                                                  #
# SEQUENCE EXTRACTION FROM SUBJECT AND QUERY FILES #
#       WITHOUT STOP CODONS WITHIN SEQUECE         #
#                                                  #
#      INPUT - TWO FILES WITH DNA SEQUENCES        # 
#              (QUERY AND SUBJECT)                 #
#                       AND                        #
#    PARSED BLAST REPORT BY tcl_blast_parser_123   #
#                                                  #
#     NOTE! BLAST REPORT CORRESPONDS TO BLASTX     #
#    (TRANSLATED DNA VERSUS PROTEIN SEQUENCES)     #
# IT MEANS THAT SUBJECT DNA FILE SHOULD CORRESPOND #
# EXACTLY TO PROTEIN QUERY FILE FOR PROPER RESULTS #
#                                                  #
####################################################

####################################################
#                                                  #
#      COPYRIGHT, ALEXANDER KOZIK, 2003 - 2007     #
#                                                  #
####################################################

def Seqs_Nficator(sequence):


	non_atgc_list = [ 'B', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'M', \
			'O', 'P', 'Q', 'R', 'S', 'U', 'V', 'W', 'X', 'Y', 'Z' ]

	have_seqs = string.upper(sequence)
	for dummy_letter in non_atgc_list:
		sequence = re.sub(dummy_letter, "N", sequence)

	return sequence

def Seqs_Extractor(in_name1, in_name3, out_name):

	global stop_codon_status

	print "INPUT FILE 1 (QR_SEQS) :   " + in_name1
	print "INPUT FILE 3 (BLASTX)  :   " + in_name3
	print "OUTPUT FILE (SEQ FRAGM):   " + out_name

	seq5_name = out_name + '.x_seq5'
	seq3_name = out_name + '.x_seq3'
	seqX_name = out_name + '.x_seqSummary'
	seqY_name = out_name + '.x_seqCDS_UP_Clean'
	seqZ_name = out_name + '.x_seqCDS_UP_All'
	seqN_name = out_name + '.x_no_hits_found'

	in_file1  = open(in_name1, "rb")
	in_file3  = open(in_name3, "rb")
	out_file  = open(out_name + '.x_seqCDS_FR_Region', "wb")
	seq5_file = open(seq5_name, "wb")
	seq3_file = open(seq3_name, "wb")
	seqX_file = open(seqX_name, "wb")
	seqY_file = open(seqY_name, "wb")
	seqZ_file = open(seqZ_name, "wb")
	seqN_file = open(seqN_name, "wb")

	### HEADER IN INFO FILE ###
	seqX_file.write('QUERY_ID' + '\t' + 'SUBJECT_ID' + '\t' + 'QR_LEN' + '\t' + 'SB_LEN' + '\t' + \
			'***' + '\t' + 'S_START' + '\t' + 'S_END' + '\t' + '***' + '\t' 'QR_ALN' + '\t' \
			+ 'SB_ALN' + '\t' + '***' + '\t' + 'Q_5_DED' + '\t' 'Q_5_REG' + '\t' + 'Q_START' + \
			'\t' + 'LENGTH' + '\t' + 'Q_END' + '\t' + 'Q_3_REG' + '\t' + 'Q_3_DED' + '\t' + 'DIRECT' + '\t' + 'CASE' + '\n')

	seq_list  = []
	hit_list  = []
	seq_array = {}
	hit_array = {}

	##########################################
	# READ SEQUENCE DATA INTO ARRAY : STEP 1 #
	#         QUERY DNA FASTA FILE           #
	##########################################
	line_counter = 0
	#################################
	while 1:
		### LAST SEQUENCE IN FASTA FILE ###
		t = in_file1.readline()
		if t == '':
			###  SUB_SEQ FUNCTION  ###
			have_seqs = "".join(my_seqs)
			seqs_len = len(have_seqs)
			seq_array[proper_id] = have_seqs
			####  END OF SUB_SEQ  ####
			break
		### SEQ PROCESSING ###
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]

		fasta_match = t[0:1]
		if fasta_match == ">":
			gi_test = t[0:4]
			if gi_test == ">gi|":
				# print gi_test
				descr_line = t
				descr_line = re.sub("^>gi\|", "", descr_line)
				descr_line = re.sub("\|", '\t', descr_line, 1)
				# print line_counter
				line_counter += 1
			else:
				descr_line = t
				descr_line = re.sub("^>", "", descr_line)
				descr_line = re.sub(" ", '\t', descr_line, 1)
				# print line_counter
				line_counter += 1
			good_head = string.split(descr_line, '\t')[0]
			try:
				long_tail = string.split(descr_line, '\t')[1]
			except:
				long_tail = "no description found"
			if good_head in seq_list:
				running_text = "\
\n\n  Ooops... ID duplication  \n\n  check input for duplications  \n\n\n ID: " + good_head + "\n\n\n"
				print running_text
				###
				### INSERT TKINTER TEXT MESSAGE BOX
				###
				break
			seq_list.append(good_head)
			print `line_counter` + '\t' + good_head
			if line_counter != 1:
				###  SUB_SEQ FUNCTION  ###
				have_seqs = "".join(my_seqs)
				seqs_len = len(have_seqs)
				seq_array[proper_id] = have_seqs
				#print good_head
				#print have_seqs
				####  END OF SUB_SEQ  ####
			have_seqs = ""
			my_seqs = []
		if fasta_match != ">" and fasta_match != "":
			proper_id = good_head
			good_name = long_tail
			my_seqs.append(t)
	#################################
	###       END OF STEP 1         #
	#################################

	#################################
	hit_count = 0
	no_hits_count = 0
	#################################
	#       READ BLASTX DATA        #
	#################################
	while 1:
		t = in_file3.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		t = t.split('\t')
		#################
		id  = t[0]
		# print 'id' + '\t' + id
		hit = t[1]
		# print 'hit' + '\t' + hit
		ds  = t[2]
		# print 'ds' + '\t' + ds
		if ds == 'no_hits_found':
			no_hits_count = no_hits_count + 1
			print "no hits found:  " + `no_hits_count`
			sequence = seq_array[id]
			sequence = Seqs_Nficator(sequence)
			find_crap = sequence.rfind('N')
			if find_crap == -1:
				crap_status = "CLEAN__NO_NN"
			if find_crap != -1:
				crap_status = "DIRTY_WITH_N"
			seqN_file.write('>' + id + '  _no_hits_found_ ' + crap_status + '\n')
			seqN_file.write(sequence + '\n')
		if ds != 'no_hits_found':
			aln_idn = t[4]
			hnb = t[7]
			hnb = int(hnb)
			# print 'hnb' + '\t' + hnb
			fr  = t[8]
			# print 'fr' + '\t' + fr
			qst = t[9]
			qst = int(qst)
			# print 'qst' + '\t' + qst
			qen = t[10]
			qen = int(qen)
			# print 'qen' + '\t' + qen
			hst = t[11]
			hst = int(hst)
			# print 'hst' + '\t' + hst
			hen = t[12]
			hen = int(hen)
			# print 'hen' + '\t' + hen
			qsl = t[13]
			qsl = qsl.split('/')
			q_len = qsl[0]
			s_len = qsl[1]
			q_len = int(q_len)
			s_len = int(s_len)
			gaps = t[14]
			#################################
			#    SUBSEQUENCE EXTRACTION     #
			#################################
			if hnb == 1:
				hit_count += 1
				print `hit_count` + '\t' + id
				########################
				sequence = seq_array[id]
				sequence = Seqs_Nficator(sequence)
				# print sequence
				seq_len  = len(sequence)
				if seq_len != q_len:
					print "__ DATA MISMATCH __"
					print "CHECK SEQUENCE: " + id
					sys.exit()
				if qen > qst:
					direction = 'forward'
					sub_seq  = sequence[(qst-1):(qen-0)]
					sub_len  = len(sub_seq)
					sub_seq5 = sequence[:(qst-1)]
					sub_seq3 = sequence[(qen-0):]
				if qen < qst:
					direction = 'reverse'
					t_r = list(sequence)	# string -> list of chars
					t_r.reverse()	# inplace reverse the list
					t_r = string.join(t_r, '')	# list of strings -> string
					t_r = string.upper(t_r)		# all uppercase
					t_r = re.sub("A", "t", t_r)	# A -> t
					t_r = re.sub("T", "a", t_r)	# T -> a
					t_r = re.sub("G", "c", t_r)	# G -> c
					t_r = re.sub("C", "g", t_r)	# C -> g
					t_r = string.upper(t_r)		# back to uppercase
					sequence = t_r
					sub_seq  = sequence[(seq_len-qst):(seq_len-qen + 1)]
					sub_len  = len(sub_seq)
					sub_seq5 = sequence[:(seq_len-qst)]
					sub_seq3 = sequence[(seq_len-qen + 1):]
				# print sub_seq
				mod_3 = math.fmod(sub_len, 3)
				## PROTEIN POSITION INTO DNA POSITION ###
				hst = hst*3 - 2
				hen = hen*3
				if mod_3 != 0:
					print "MOD IS NOT 3 !!!"
					time.sleep(1)
				n_find  = sub_seq.rfind('N')
				################################
				acase = "WHATEVER"
				seq_len2 = s_len*3
				sub_len2 = hen - hst + 1
				################################
				if direction == 'forward':
					qs_real = qst
					q5_reg  = len(sub_seq5)
					q5_ded  = qst - hst
					qe_real = qen
					q3_reg  = len(sub_seq3)
					q3_ded  = (seq_len - qen) - (seq_len2 - hen)
				if direction == 'reverse':
					qs_real = seq_len-qst
					q5_reg  = len(sub_seq5)
					q5_ded  = (seq_len-qst) - hst
					qe_real = seq_len-qen + 1
					q3_reg  = len(sub_seq3)
					q3_ded  = (seq_len - (seq_len-qen + 1)) - (seq_len2 - hen)
				################################
				if q5_ded >= 0 and q3_ded >= 0:
					acase = '4'
				if q5_ded < 0 and q3_ded < 0:
					acase = '1'
				if q5_ded >= 0 and q3_ded < 0:
					acase = '2'
				if q5_ded < 0 and q3_ded >= 0:
					acase = '3'
				################################
				if n_find == -1:
				################################
					out_file.write('>' + id + '.CDS' + ' length: ' + `sub_len` + ' ' + 'mod: ' + `mod_3` + ' ' + \
					direction + ' ' + `qst` + '-' + `qen` + \
					' [' + hit + ' ' + `hst` + '-' + `hen` + ']' + '\n' + sub_seq + '\n')
				################################
					if len(sub_seq3) >= 12:
						seq3_file.write('>' + id + '.3PR' + ' 3-REGION ' + direction + '\n' + sub_seq3 + '\n')
					if len(sub_seq5) >= 12:
						seq5_file.write('>' + id + '.5PR' + ' 5-REGION ' + direction + '\n' + sub_seq5 + '\n')
				################################
					seqX_file.write(id + '\t' + hit + '\t' + `seq_len` + '\t' + `seq_len2` + '\t' + '***' + '\t' + \
						`hst` + '\t' + `hen` + '\t' + '***' + '\t' + \
						`sub_len` + '\t' + `sub_len2` + '\t' + '***' + '\t' + \
						`q5_ded` + '\t' + `q5_reg` + '\t' + `qs_real` + '\t' + `seq_len` + '\t' + \
						`qe_real` + '\t' + `q3_reg` + '\t' + `q3_ded` + '\t' + direction + '\t' + acase + '\n')
				################################
					sub_seq3_low = string.lower(sub_seq3)
					sub_seq5_low = string.lower(sub_seq5)
					seqY_file.write('>' + id + ' ' + direction + ' ' + `qst` + '-' + `qen` + \
								' [' + hit + ' ' + `hst` + '-' + `hen` + ']' + \
								'\n' + sub_seq5_low + sub_seq + sub_seq3_low + '\n')
				################################
				sub_seq3_low = string.lower(sub_seq3)
				sub_seq5_low = string.lower(sub_seq5)
				seqZ_file.write('>' + id + ' ' + direction + ' ' + `qst` + '-' + `qen` + \
							' [' + hit + ' ' + `hst` + '-' + `hen` + ']' + \
							'\n' + sub_seq5_low + sub_seq + sub_seq3_low + '\n')
				################################

	in_file1.close()
	in_file3.close()
	out_file.close()
	seq5_file.close()
	seq3_file.close()
	seqX_file.close()
	seqY_file.close()
	seqZ_file.close()
	seqN_file.close()

##################
#                #
#   MAIN BODY    #
#                #
##################

import math
import re
import sys
import string
import time
import os

if __name__ == "__main__":
	if len(sys.argv) <= 3 or len(sys.argv) > 4:
		print "Program usage: "
		print "input_file1(query_seqs) input_file3(blastx.info2) output_file"
		sys.exit()
	if len(sys.argv) == 4:
		in_name1  = sys.argv[1]
		# in_name2  = sys.argv[2]
		in_name3  = sys.argv[2]
		out_name  = sys.argv[3]
		# out_dir   = sys.argv[5]
		Seqs_Extractor(in_name1, in_name3, out_name)
### THE END ###
