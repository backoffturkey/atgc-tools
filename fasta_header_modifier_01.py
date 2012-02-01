#!/usr/bin/python

########################
# COPYRIGHT 2010       #
# Alexander Kozik      #
# http://www.atgc.org/ #
# akozik@atgc.org      #
########################

def Fasta_Modifier(in_name, out_name1, out_name2, group_file, tab_format):

	print in_name + ' ' + out_name1 + ' ' + out_name2

	in_file  = open(in_name,  "rb")
	out_file1 = open(out_name1, "wb")
	out_file2 = open(out_name2, "wb")

	## ANNOTATION ARRAY ##
	frame_array = {}
	
	mod_val = 1000
	# mod_val = 100000
	id_count = 0
	all_count = 0
	extr_count = 0
	except_count = 0
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
			id = u[0]
			
			if tab_format == "TAB-3":
				hit_id = u[1]
				description = u[2]
				frame_array[id] = " {HIT_ID: " + hit_id + " DESCR: "  + description + "} " 
			
			if tab_format == "TAB-CLC-E":
				n_hits   = u[1]
				e_value  = u[2]
				e_hit_id = u[3]
				e_descr  = u[4]
				
				description = " {HIT_ID: " + e_hit_id + " EXP:(" + e_value + ")" + " DESCR: " + e_descr + "} " 
				
				if n_hits == "0":
					description = " {_no_hits_found_} "
				## ARRAY WITH ANNOTATION ##
				frame_array[id] = description
			
			## ID COUNT ##
			id_count = id_count+1
			print_update = math.fmod(id_count, mod_val)
			if print_update == 0:
				print `id_count` + " IDs in group file processed "
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
					annotation = frame_array[proper_id]
					sub_seqs = have_seqs
					out_file1.write('>' + proper_id + ' ' + good_name + annotation + '\n')
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
						annotation = frame_array[proper_id]
						sub_seqs = have_seqs
						out_file1.write('>' + proper_id + ' ' + good_name + annotation + '\n')
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
		print "[input_file] [output_file1] [output_file2] [id_description] [tab_format]"
		print ""
		print "Sequences listed in id_description table will be extracted and annotated into FASTA output_file1"
		print "Sequences NOT listed in id_description table will be written into output_file2"
		print ""
		sys.exit()
	if len(sys.argv) == 6:
		in_name    = sys.argv[1]
		out_name1  = sys.argv[2]
		out_name2  = sys.argv[3]
		group_file = sys.argv[4]
		tab_format = sys.argv[5]
		if tab_format != "TAB-3" and tab_format != "TAB-CLC-E":
			print "tab_format must be \'TAB-3\' or \'TAB-CLC-E\'"
			sys.exit()
		if in_name != out_name1 and in_name != out_name2:
			Fasta_Modifier(in_name, out_name1, out_name2, group_file, tab_format)
		else:
			print "Output should have different name than Input"
			sys.exit()
