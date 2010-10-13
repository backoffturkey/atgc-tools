#!/usr/bin/expect -f

#####################################
#                                   #
#           MAIN BODY               #
#                                   #
#####################################

# standard input/output dialog #

puts -nonewline "Enter the SOURCE (GenBank) file name: "
flush stdout

set sourcefile [gets stdin]

while {[catch {set f_in [open "$sourcefile"  "r"]}] == 1} {
        puts stderr "Error opening $sourcefile"
        puts -nonewline "Please re-enter the SOURCE file name: "
        flush stdout
        set sourcefile [gets stdin]
}

puts -nonewline "Enter the DESTINATION (Fasta) file name: "
flush stdout

set destfile [gets stdin]

set f_out [open "$destfile"  "w"]
set f_dic [open "$destfile\.dict" "w"]

# end of input/output dialog #

# PATTERNS #

set lineNumber 0
set first_line_to_process 0

set locus_p      "LOCUS       "
set definition_p "DEFINITION  "
set accession_p  "ACCESSION   "
set version_p    "VERSION     "
set organism_p   "  ORGANISM  "
set authors_p    "  AUTHORS   "
set origin_p     "ORIGIN      "
set end_data_p   "//"
set strain_p     "                     /strain="
set cultivar_p   "                     /cultivar="
set ecotype_p    "                     /ecotype="

set locus_s      "WHATEVER"
set definition_s "WHATEVER"
set accession_s  "WHATEVER"
set version_s    "WHATEVER"
set organism_s   "WHATEVER"
set origin_s     "WHATEVER"
set seqs_s       "WHATEVER"
set end_data_s   "WHATEVER"

set locus_c 0
set good_c  0

while {[gets $f_in  current_line] >= 0}  {

        if {$lineNumber >= $first_line_to_process} {

	    set t [string range $current_line 0 11]
	    set u [string range $current_line 0 1 ]
	    set x $current_line

	    set t_s_28 [string range $current_line 0 28]
	    set t_c_30 [string range $current_line 0 30]
	    set t_e_29 [string range $current_line 0 29]

	    if { $t == $locus_p } {
		regsub {^LOCUS       } $x "" locus_f 
		regsub { .*} $locus_f "" locus_f
		set locus_s      "DONE"
		set definition_s "WHATEVER"
		set accession_s  "WHATEVER"
		set version_s    "WHATEVER"
		set organism_s   "WHATEVER"
		set origin_s     "WHATEVER"
		set seqs_s       "WHATEVER"
		set end_data_s   "WHATEVER"
		incr locus_c
		# puts "$locus_c    $locus_f"
		### UNSET ALL VARIABLES ###
		set short_name   "SHORT_NAME_TEST"
		set accession_f  "ACCESSION_TEST"
		set version_f    "VERSION_TEST"
		set genbank_id   "GENBANK_ID_TEST"
		set definition_f "DEFINITION_TEST"
		set organism_f   "ORGANISM_TEST"
		set authors_f    "AUTHORS_TEST"
		set strain_f     "-=ND4=-"
		set cultivar_f   "-=ND5=-"
		set ecotype_f    "-=ND6=-"
	    }

	    if { $t == $definition_p } {
		regsub {^DEFINITION  } $x "" definition_f
		regsub { \[.*} $definition_f "" definition_f
		regsub -all {  } $definition_f " " definition_f
		set definition_s "DONE"
		# puts "$locus_c | $locus_f | $definition_f"
	    }

	    if { $t == $accession_p } {
		regsub {^ACCESSION   } $x "" accession_f
		regsub { .*} $accession_f "" accession_f
		set accession_s "DONE"
	    }

	    if { $t == $version_p } {
		regsub {^VERSION     } $x "" version_f
		regsub {.* GI\:} $x "" genbank_id
		regsub { .*} $version_f "" version_f
		set version_s "DONE"
	    }

	    if { $t == $organism_p } {
		regsub {^  ORGANISM  } $x "" organism_f
		regsub -all {\'} $organism_f "" organism_f
		set organism_a [split $organism_f " "]
		set name1 [lindex $organism_a 0]
		set name2 [lindex $organism_a 1]
                if { [llength $organism_a] < 2 } {  
                    set name2 "XXX"
                    puts "NAME IS SHORT!"
                }

		set name1s [string range $name1 0 3]
		set name2s [string range $name2 0 3]
		set short_name "$name1s\_$name2s"
		# regsub {\'} $short_name "" short_name
		regsub -all {\.} $short_name "" short_name
		regsub -all {\"} $short_name "" short_name
		regsub -all {\:} $short_name "" short_name
		set organism_s "DONE"
	    }

            if { $t == $authors_p } {
                regsub {^  AUTHORS   } $x "" authors_f
                }

            #### STRAIN ####
            if { $t_s_28 == $strain_p } {
		regsub {^                     /strain=} $x "" strain_f
               }
            if { $t_c_30 == $cultivar_p } {
		regsub {^                     /cultivar=} $x "" cultivar_f
               }
            if { $t_e_29 == $ecotype_p } {
		regsub {^                     /ecotype=} $x "" ecotype_f
               }

	    if { $t == $origin_p } {
		set seqs_s "TPAXHEM"
		# set name_length [ expr [string length $short_name] + [string length $locus_f] + 1 ]
		set name_length [ expr [string length $short_name] + [string length $accession_f] + 1 ]
		set warning_mess ""
		if { $name_length > 24 } {
			set warning_mess " _LONG_NAME_WARNING_ "
			puts " LONG NAME:  $short_name\.$locus_f"
			after 2000
		}
		puts "$locus_c | $locus_f | $accession_f | $version_f | $genbank_id | $organism_f | $short_name | NAME_LEN: $name_length"
		puts $f_out "\>$short_name\.$accession_f$warning_mess \( $version_f GI\:$genbank_id \) \{ $definition_f \} \[ $organism_f \]"
		puts $f_dic "$short_name\t |1| \t$organism_f\t |2| \t$short_name\.$accession_f\t |3| \t$definition_f\t |4| \t$strain_f\t |5| \t$cultivar_f\t |6| \t$ecotype_f\t |7| \t -=| $authors_f |=- "
		# puts $f_out "\>$short_name\.$locus_f$warning_mess \( $version_f GI\:$genbank_id \) \{ $definition_f \} \[ $organism_f \]"
		# puts $f_dic "$short_name\t |1| \t$organism_f\t |2| \t$short_name\.$locus_f\t |3| \t$definition_f\t |4| \t$strain_f\t |5| \t$cultivar_f\t |6| \t$ecotype_f\t |7| \t -=| $authors_f |=- "
	    }

	    if  { $u == $end_data_p } {
		set seqs_s "XPEH BAM"
		puts -nonewline $f_out "\n"
		set locus_s "WHATEVER"
	    }

	    if { $seqs_s == "TPAXHEM" && $t != $origin_p && $locus_s == "DONE" && $definition_s == "DONE" && $accession_s == "DONE" && $version_s == "DONE" && $organism_s == "DONE" } {
		regsub -all { } $x "" dummy_seqs
		regsub -all {[0-9]} $dummy_seqs "" dummy_seqs
		set dummy_seqs [string toupper $dummy_seqs]
		# puts $dummy_seqs
		puts -nonewline $f_out $dummy_seqs
	    } 
        }
        incr lineNumber
	# puts $lineNumber
    }

close $f_in
close $f_out
close $f_dic

