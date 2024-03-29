#!/usr/bin/python

########################
# COPYRIGHT 2010       #
# Alexander Kozik      #
# http://www.atgc.org/ #
# akozik@atgc.org      #
########################

def Seqs_Drobilka(in_name, out_name1, out_name2, group_file, seqs_range):

	print in_name + ' ' + out_name1 + ' ' + out_name2

	in_file  = open(in_name,  "rb")
	out_file1 = open(out_name1, "wb")
	out_file2 = open(out_name2, "wb")

	frame_array = {}
	start_array = {}
	end_array   = {}

	mod_val = 100000
	id_count = 0
	all_count = 0
	extr_count = 0
	except_count = 0
	print_frame = "FALSE"
	### TRY TO READ GROUP LIST ###
	try:
		frame_file = open(group_file, "rb")
		print "USING GROUP LIST"
		while 1:
			u = frame_file.readline()
			if u == '':
				break
			if '\n' in u:
				u = u[:-1]
			if '\r' in u:
				u = u[:-1]
			u = u.split('\t')
			fl = u[0]
			frame_array[fl] = fl
			id_count = id_count+1
			print_update = math.fmod(id_count, mod_val)
			if print_update == 0:
				print `id_count` + " IDs in group file processed "
			if seqs_range == "RANGE":
				start_array[fl] = int(u[1])
				end_array[fl]   = int(u[2])
	except:
		print "DID NOT FIND GROUP FILE:  " + group_file
		sys.exit()

	print "------------------------------------------------"
	print `id_count` + " IDs in group file were found "
	print "------------------------------------------------"
	time.sleep(2)
	
	fasta_id_array = []
	line_counter = 0
	have_seqs = ""
	proper_id = ""
	my_seqs = []

	while 1:
		t = in_file.readline()
		if t == '':
			###  SUB_SEQ FUNCTION  ###
			have_seqs = "".join(my_seqs)
			seqs_len = len(have_seqs)
			###   STRING PROCESSING   ###
			if seqs_len != 0:
				have_seqs = "".join(my_seqs)
				seqs_len = len(have_seqs)

			if have_seqs != "":
				try:
					frame_array[proper_id]
					if seqs_range == "RANGE":
						sub_seqs = have_seqs[start_array[proper_id]:end_array[proper_id]]
						out_file1.write('>' + proper_id + ' ' + good_name + ' FRAGMENT ' + `start_array[proper_id]` + ':' + `end_array[proper_id]` + '\n')
					if seqs_range == "FULL":
						sub_seqs = have_seqs
						out_file1.write('>' + proper_id + ' ' + good_name + '\n')
					out_file1.write(sub_seqs + '\n' + '\n')
					extr_count = extr_count + 1
					print_update = math.fmod(extr_count, mod_val)
					if print_update == 0:
						print `extr_count` + " IDs were extracted ... "
				except:
					### print proper_id  " NOT EXTRACTED "
					# mooba = "booba"
					except_count = except_count + 1
					out_file2.write('>' + proper_id + ' ' + good_name + '\n')
					out_file2.write(have_seqs + '\n' + '\n')
				break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]

		fasta_match = t[0:1]
		if fasta_match == ">":
			descr_line = t
			descr_line = re.sub("^>", "", descr_line)
			descr_line = re.sub(" ", '\t', descr_line, 1)
			good_head = string.split(descr_line, '\t')[0]
			line_counter += 1
			try:
				long_tail = string.split(descr_line, '\t')[1]
			except:
				long_tail = ""
			if line_counter != 1:
				have_seqs = "".join(my_seqs)
				seqs_len = len(have_seqs)

				if have_seqs != "":
					try:
						frame_array[proper_id]
						if seqs_range == "RANGE":
							sub_seqs = have_seqs[start_array[proper_id]:end_array[proper_id]]
							out_file1.write('>' + proper_id + ' ' + good_name + ' FRAGMENT ' + `start_array[proper_id]` + ':' + `end_array[proper_id]` + '\n')
						if seqs_range == "FULL":
							sub_seqs = have_seqs
							out_file1.write('>' + proper_id + ' ' + good_name + '\n')
						out_file1.write(sub_seqs + '\n' + '\n')
						extr_count = extr_count + 1
						print_update = math.fmod(extr_count, mod_val)
						if print_update == 0:
							print `extr_count` + " IDs were extracted ... "
					except:
						### print proper_id  " NOT EXTRACTED "
						# mooba = "booba"
						except_count = except_count + 1
						out_file2.write('>' + proper_id + ' ' + good_name + '\n')
						out_file2.write(have_seqs + '\n' + '\n')

			have_seqs = ""
			my_seqs = []
			
		if fasta_match != ">" and fasta_match != "":
			proper_id = good_head
			good_name = long_tail
			my_seqs.append(t)

	in_file.close()
	out_file1.close()
	out_file2.close()

	print ""
	print " PROCESSING DONE! "
	print `extr_count` + " - total number of extracted IDs "
	print `except_count` + " - total number of exceptions "
	print `line_counter` + " - all items"
	print ""

import time
import math
import re
import sys
import string
if __name__ == "__main__":
	if len(sys.argv) <= 5 or len(sys.argv) > 6:
		print ""
		print "Program usage: "
		print "[input_file] [output_file1] [output_file2] [group_id_file] [seqs_range]"
		print ""
		print "Sequences listed in group_file will be extracted according to FASTA header IDs into output_file1"
		print "Sequences NOT listed in group_file will be written into output_file2"
		print ""
		exit
	if len(sys.argv) == 6:
		in_name    = sys.argv[1]
		out_name1  = sys.argv[2]
		out_name2  = sys.argv[3]
		group_file = sys.argv[4]
		seqs_range = sys.argv[5]
		if seqs_range != "FULL" and seqs_range != "RANGE":
			print "seqs_range must be \'FULL\' or \'RANGE\'"
			sys.exit()
		if in_name != out_name1 and in_name != out_name2:
			Seqs_Drobilka(in_name, out_name1, out_name2, group_file, seqs_range)
		else:
			print "Output should have different name than Input"
			sys.exit()
