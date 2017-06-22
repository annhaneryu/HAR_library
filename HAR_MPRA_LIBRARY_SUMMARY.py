# SUMMARY OF HAR MPRA LIBRARY DESIGN #

#################################################################################################################
# MODIFY hars_merged_hg19.bed file for extracting MAF and FASTA files and generating permutations
#################################################################################################################

# modify the hars_merged_hg19.bed file such that the HARs have a new start site, one less than the original start site
# keep the 2xHAR coordinates the same

directory = '/Users/haneryu/Documents/Pollard_Lab_2013/'

file1 = open(directory + 'hars_merged_hg19.bed', 'rb')
file2 = open(directory + 'zero_based_HARs_hg19.bed', 'wb')
file3 = open(directory + 'zero_based_2xHARs_hg19.bed', 'wb')

HAR_dictionary = {}
TWOxHAR_dictionary = {}

HARs = []
TWOxHARs = []
HAR_lines = []
TWOxHAR_lines = []

for line in file1:
	#print line
	TWOxHAR_lines.append(line)
	data = line.klklip('\r\n').split('\t')
	chrom = data[0]
	start = data[1]
	end = data[2]
	name = data[3]
	if '2x' not in name:
		print name
		HARs.append(name)
		new_start = int(start) - 1 
		print start
		print new_start
		new_line = [chrom, new_start, end, name]
		print new_line
		HAR_lines.append(new_line)
		HAR_dictionary[name] = new_line
	else:
		#print name
		TWOxHARs.append(name)
		TWOxHAR_dictionary[name] = line

for HAR in HARs:
	new_line = HAR_dictionary[HAR]
	print new_line
	#file2.write(new_line + '\n')
	file2.write('\t'.join(map(str, new_line)) + '\n')
file2.close()

for TWOxHAR in TWOxHARs:
	line = TWOxHAR_dictionary[TWOxHAR]
	print line
	file3.write(line)
file3.close()


#################################################################################################################
# FIRST STEP. INPUT BED FILE DIRECTORY --> EXTRACTED MAF FILE PER HAR FROM BED FILES
#################################################################################################################

# BEFORE EXTRACTING FASTAs per HAR
#to extract MAF files per HAR BED file
import sys
import os
import re
import subprocess

#dir='/home/hane/HAR_BEDfiles/'
dir = '/home/hane/HAR_BEDfiles/All_MUTANT_OFF_BY_ONE_BEDS'

#to add one to the start site of all BED files in this BED file directory, use this command:
# ~/HAR_BEDfiles/All_BEDfiles$ for fff in *.bed ; do echo "PROCESSING $fff" ; cat $fff | awk '{print $1 "\t" $2-1 "\t" $3 "\t"$4}' > ../All_MUTANT_OFF_BY_ONE_BEDS/mutant.$fff ; done

for bed in os.listdir(dir):
	if bed.endswith('.bed'):
		print bed
		bed_name = bed.split('_')
		chrom_name = bed_name[1]
		#chrom_name = re.match('chr\d+', bed).group(0)
#		if chrX, then use a break statement
		maf_file = chrom_name + '.maf'
		output=bed.replace('.bed', '.maf').replace(';', '\;')
		print "hane is awesome"
		command="/db/projects/Common/bin/mafsInRegion " + dir + '/' + bed + ' ' + '/home/hane/HAR_BEDfiles/MAF_per_HAR' + '/'  + output + ' ' + '/db/projects/hane/ucsc_multiz46_maf/MAF_per_chrom'  + '/' + maf_file
		command = command.replace(';', '\;')
		print command
		print output
		result=subprocess.call(command, shell=True)
		if result!=0:
			sys.exit()

#################################################################################################################
# SECOND STEP. INPUT MAF FILE DIRECTORY --> EXTRACTED FASTA FILE PER HAR FROM MAF FILES
#################################################################################################################

#to extract FASTA files per HAR MAF file using msa_view

import sys
import os
import re
import subprocess

#msa_view maf --in-format MAF --out-format FASTA --seqs hg19, panTro2 > output

dir = '/home/hane/HAR_MPRA_library_design/MAF_per_HAR'
out_dir = '/home/hane/HAR_MPRA_library_design/FASTA_per_HAR/'

for maf in os.listdir(dir):
	if maf.endswith('.maf'):
		print maf
		output = maf.replace('.maf', '.fa')
		#command = msa_view maf --in-format MAF --out-format FASTA --seqs hg19, panTro2 > out_dir + output
		command = '/home/Common/opt/Bio/phast/bin/msa_view' + ' ' +  dir + '/' + maf + ' ' + '--in-format MAF --seqs hg19,panTro2 > ' + out_dir + output
		print command
		print output
		result = subprocess.call(command, shell=True) 
		print "yes it worked!"
		if result!=0:
			sys.exit()


#################################################################################################################
# THIRD STEP. REMOVE GAP CHARACTERS IN FASTA FILES EXTRACTED FROM MAF FILES
#################################################################################################################

import os

def load_sequences(fasta_file):
		file1 = open(fasta_file, 'rb')
		data = file1.readlines()
		str_data = ''.join(data)
		seq_data = filter(None, str_data.split('>'))

		hum_data = filter(None, seq_data[0].split('\n'))
		hum_seq = ''.join([x.strip() for x in hum_data[1:]])

		chimp_data = filter(None, seq_data[1].split('\n'))
		chimp_seq = ''.join([x.strip() for x in chimp_data[1:]])

		return hum_seq, chimp_seq

def remove_gaps(hum_seq, chimp_seq, out_file):
		gap_positions = []
		new_human_seq = []
		new_chimp_seq = []
		human = '> hg19'
		chimp = '> panTro2'
		for i in range(0, len(hum_seq)):
			chimp_char = chimp_seq[i]
			hum_char = hum_seq[i]
			if (chimp_char != '-') and (hum_char != '-'):
				new_human_seq.append(hum_char)
				new_chimp_seq.append(chimp_char)
				#print chimp_char, hum_char, i
				#if chimp_char = hum_char:
				#	gap_position.append(i)
			else:
				gap_positions.append(i)
		hg19_seq = ''.join(new_human_seq)
		panTro2_seq = ''.join(new_chimp_seq)
		out_file.write('\n'.join([human, hg19_seq, chimp, panTro2_seq]))
		out_file.close()

# def remove_gaps(hum_seq, chimp_seq, gap_positions, out_file):
# 		for i in range(0, len(hum_seq)):
# 			chimp_char = chimp_seq[i]
# 			hum_char = hum_seq[i]

directory = '/home/hane/HAR_FASTA_files/modified_FASTA_per_HAR'
output_directory = '/home/hane/HAR_FASTA_files/modified_FASTA_per_HAR/further_modified_HAR_FASTAs'
fasta_files = os.listdir(directory)
for fasta_file in fasta_files:
		print fasta_file
		if fasta_file.endswith('.fa'):
				hum_seq, chimp_seq = load_sequences(directory + '/' + fasta_file)
				out_file_name = ('new' + '_' + fasta_file)
				out_file = open(output_directory + '/' + out_file_name, 'wb')
				remove_gaps(hum_seq, chimp_seq, out_file)              
#				out_file.write('\t'.join(map(str, [hum_char, chimp_char, i])) + '\n')	


#################################################################################################################
# FOURTH STEP. SPLIT UP FASTA FILE PER HAR 
#################################################################################################################

#!/usr/bin/python

import pdb # debugger : run the script like so: python -m pdb hane_find_forbidden_sequences.py 
 # type "c" to continue execution (you have to do this when you start or the program won't run... annoying

import os
import subprocess
import sys
import re
from itertools import chain

input_fasta_dir = "/home/hane/MODIFIED_FASTA_files"
#input_fasta_dir = 'whatever directory has all 2xHAR- and modified HAR- FASTA files'
output_directory = "/home/hane/MODIFIED_SPLIT_UP_HARs"

if not os.path.exists(output_directory):
	os.makedirs(output_directory)
	pass

def load_sequences(fasta_file):
	file1 = open(fasta_file, 'rb')
	data = file1.readlines()
	#
	str_data = ''.join(data)
	seq_data = filter(None, str_data.split('>'))
	#
	hum_data = filter(None, seq_data[0].split('\n'))
	hum_seq = ''.join([x.strip() for x in hum_data[1:]])
	#
	return hum_seq
#
#
MAX_SEQ_LEN = 170
#too_long = "\w{170}" # any "word" character 170 times
XbaI_F = 'TCTAGA'
XbaI_R = 'AGATCT'
EcoRI_F = 'GAATTC'
EcoRI_R = 'CTTAAG'
SbfI_F = 'CCTGCAGG'
SbfI_R = 'GGACGTCC'
Four_Gs = 'GGGG'  # this can still generate 4 gs in a row in unusual circumstances, but it doesn't really matter since it requires a ton of Gs in a row and will probably never happen

lists = [EcoRI_F, EcoRI_R, SbfI_F, SbfI_R, XbaI_F, XbaI_R, Four_Gs]

good_HARS = []

fasta_files = os.listdir(input_fasta_dir)
for fasta_file in fasta_files:
	#
	if (not fasta_file.endswith('.fa')):
		print "SKIPPING THE APPARENTLY NON-FASTA FILE" + fasta_file + "..." + " ----------------------------------------------------------- "
		continue # next iteration of the loop after skipping this weird file
	#
	#
	#fasta_file_name = fasta_file.strip('.fa').split('_')
	#HAR_name = fasta_file_name[2]
	control_name = fasta_file.strip('.fa')
	hum_seq = load_sequences( os.path.join(input_fasta_dir, fasta_file) )
	har_problems = []
	har_problem_starts = []
	splitLocationsWithDuplicatesArray = []
        cutsBeforeHandlingLengthSet    = []
	Four_Gs_sites = []
	new_HAR_hg_seqs = []
	#
	for illegalSeq in lists:
		if illegalSeq in hum_seq:
			#print fasta_file, illegalSeq
			har_problems.append(illegalSeq)
			hg_Four_Gs_starts = [match.start() for match in re.finditer(Four_Gs, hum_seq)]
			Four_Gs_sites.append(hg_Four_Gs_starts)
			#hg_too_long_starts = [match.start() for match in re.finditer(too_long, hum_seq)]
			#if (len(hg_too_long_starts) > 0):
			#		print "This many human seqs in " + HAR_name + " were too long, and were cut early:", len(hg_too_long_starts)
			hg_XbaI_F_starts = [match.start() for match in re.finditer(XbaI_F, hum_seq)]
			hg_XbaI_R_starts = [match.start() for match in re.finditer(XbaI_R, hum_seq)]
			hg_EcoRI_F_starts = [match.start() for match in re.finditer(EcoRI_F, hum_seq)]
			hg_EcoRI_R_starts = [match.start() for match in re.finditer(EcoRI_R, hum_seq)]
			hg_SbfI_F_starts = [match.start() for match in re.finditer(SbfI_F, hum_seq)]
			hg_SbfI_R_starts = [match.start() for match in re.finditer(SbfI_R, hum_seq)]
			har_problem_starts.append(hg_Four_Gs_starts)
			#har_problem_starts.append(hg_too_long_starts)
			har_problem_starts.append(hg_XbaI_F_starts)
			har_problem_starts.append(hg_XbaI_R_starts)
			har_problem_starts.append(hg_EcoRI_F_starts)
			har_problem_starts.append(hg_EcoRI_R_starts)
			har_problem_starts.append(hg_SbfI_F_starts)
			har_problem_starts.append(hg_SbfI_R_starts)
			pass
		
		lotsOfStarts = list(chain(*har_problem_starts)) # somehow collapses the data structure
		splitLocationsWithDuplicatesArray = map(lambda (item): item+3, lotsOfStarts)
		cutsBeforeHandlingLengthSet = sorted(set(splitLocationsWithDuplicatesArray)) # remove duplicates!
		pass # end of the for loop
	#
	#
	cutFinalLocArr = [0] + cutsBeforeHandlingLengthSet + [len(hum_seq)] # always "cut" at the 0th base and last base (metaphorically speaking)
	#
        print cutFinalLocArr
        #
	# ===================== HERE is the part where we make sure the sequences are (after cutting illegal sequences) not any longer than MAX_SEQ_LEN
	# note that we MODIFY the cutFinalLocArr here!
	foundAnythingWeird = True
	while foundAnythingWeird:
		foundAnythingWeird = False # hopefully it's OK this time...
		for cutIndex in range(0, (len(cutFinalLocArr) - 1) ):
			cutLeft  = cutFinalLocArr[cutIndex]
			cutRight = cutFinalLocArr[cutIndex+1]
			print "Cut left: ", cutLeft, "Cut right: ", cutRight, "Size: ", (cutRight - cutLeft)
			if (cutRight - cutLeft > MAX_SEQ_LEN):
				foundAnythingWeird = True # we changed something, so we need to check again!
				#print "For HAR name " + HAR_name + ", inserted a new split at " + str(cutLeft + MAX_SEQ_LEN) + " to avoid having a sequence longer than " + str(MAX_SEQ_LEN)
				cutFinalLocArr = sorted(set(cutFinalLocArr + [cutLeft + MAX_SEQ_LEN])) # add a new element and sort the array again
				print "Updated locations array is:", cutFinalLocArr
				break # ok, now check the array from scratch again!! This could actually be a GOTO believe it or not
				pass
			pass # end "while cutIndex..." loop
		pass # end of "while foundAnythingWeird"

        print cutFinalLocArr
	# ===================== Now cutFinalLocArr should have been modified to only have sequences <= MAX_SEQ_LEN
	#
	numSequencesAfterCutting = (len(cutFinalLocArr) - 1) # Note that the 'cutFinalLocArr' has 0 and the length as 'cut' locations
	for cutIndex in range(0, (len(cutFinalLocArr) - 1) ):
		cutLeft = cutFinalLocArr[cutIndex]
		cutRight = cutFinalLocArr[cutIndex+1]
		hum_cut_seq = hum_seq[cutLeft:cutRight] # note: python slices are non-inclusive of the right side!! So for example, 0:1 just gets you the 0th element, NOT two elements!!!
		#
		cutRangeStringIndexedFromONE_not_ZERO = "range_inclusive=" + str(cutLeft+1) + "_to_" + str(cutRight) + ""
		#
		metadataStr = "cut=" + str(cutIndex+1) + "_of_" + str(numSequencesAfterCutting) + "__length=" + str(len(hum_cut_seq)) + "__" + cutRangeStringIndexedFromONE_not_ZERO
		#
		if (len(hum_cut_seq) > MAX_SEQ_LEN):
                        print "UH OH, the human sequence length didn't get properly cut, it is of length " + str(len(hum_cut_seq))
                        assert(len(hum_cut_seq) <= MAX_SEQ_LEN )
		#
		output_filename = "split_up_" + control_name + "__" + metadataStr + "__" + control_name + ".fa"
		#output_filename = "split_up_" + HAR_name + "__" + metadataStr + "__" + HAR_name + ".fa" # the har name MUST APPEAR at the very end, or the subsequent script will freak out when it can't find it. That's why the har name is in here TWO TIMES.
		output_filepath = os.path.join(output_directory, output_filename)
		out = open(output_filepath, "w")
		out.write(">" + control_name + "__" + metadataStr)
		out.write("\n")
		out.write(hum_cut_seq)
		out.write("\n")
		#print "HUMAN CUT SEQ:", hum_cut_seq
		out.close()
		#
		print "Wrote a potentially-cut sequence to the fasta file <" + output_filepath + ">"
		#
		pass # enf of the for loop

#################################################################################################################
# FIFTH STEP. CREATE THREE REQUIRED FILES FOR PERMUTATION GENERATION
#################################################################################################################

# 1) FIX THE hars_merged_hg19.bed file ----> DONE!!! use zero_based_all_HARs_hg19.bed in /home/hane/REQUIRED_FILES
# 2) recreate the SNPs_per_HARandOligos.bed file ---> DON'T NEED TO MODIFY
# 3) recreate the hane_hg19_snps_final.txt --> NEGATIVE STRAND SNP MISANNOTATION PROBLEM FIXED !!! use hane_hg19_snps_final.txt in /home/hane/REQUIRED_FILES/hane_hg19_snps_final.txt

# ALL 1067 FASTA files for SPLIT UP HARs and 2xHARs are in the same directory in lightouse: /home/hane/MODIFIED_SPLIT_UP_HARs
# ALL 284 SPLIT UP HARs are in this directory: /home/hane/MODIFIED_SPLIT_UP_HARs
# ALL SPLIT UP 783 2xHARs are in this directory: /home/hane/MODIFIED_SPLIT_UP_HARs/2xHARs

# use following script:  0b_Analyze_kp_210_hane_ryu_HAR_permutations.py in /home/hane/

#!/usr/bin/python

print "Note: this is the SECOND script that should be run, after you cut up the hars FIRST."

import pdb # debugger : run the script like so: python -m pdb 0_Analyze..........py
 # type "c" to continue execution (you have to do this when you start or the program won't run... annoying

import glob
import re
import os
import sys
import itertools
import copy
import subprocess
import sets

# How did we get the SNP_DETAILS file?
#        Answer:
#        First we downloaded "snp137.txt.gz" from the UCSC genome browser.
#        This has all the human SNP locations.
#        Then:    zcat snp137.txt.gz | grep "single" snp137.txt | cut -f 2-5,7,10 > snp137_details_only.txt
#        Now that makes a much smaller file that we can use!

# look at each file in "to process"

# get the name of the file. That tells you the SNP name

#subprocess.check_output(["ls", "-l", "/dev/null"])



# Thing 1:
# > hg19
# AACGTCCCAGAAAT
# > panTro2
# AATGTCTGAGAAAT

HAR_DIR = "/home/hane/modified_MAF_per_HAR/missing_HAR_FASTAs/modified_HAR_FASTAs/split_up_HARs" #/home/hane/MODIFIED_SPLIT_UP_HARs

OUT_DIR = "/home/hane/modified_MAF_per_HAR/missing_HAR_FASTAs/modified_HAR_FASTAs/split_up_HARs/permutations" #/home/hane/PERMUTATIONS_MODIFIED_HARS

HAR_COORDS_FILENAME = "/home/hane/REQUIRED_FILES/zero_based_all_HARs_hg19.bed" # looks like this: chr1   19688135   19688406   2xHAR.234
# gives us the coordinates for each HAR

SNP_PER_HAR_FILENAME = "/home/hane/REQUIRED_FILES/SNPs_per_HARandOligos.bed"
# lets us figure out which "rsID" is associated with each HAR
# also gives us the coordinates of the SNP!
# but: doesn't give us what the snp turns into (A/C/G/T whatever)

SNP_LETTERS_FILENAME = "/home/hane/REQUIRED_FILES/hane_hg19_snps_final.txt"

MAX_ALLOWED_DIFFERENCES = 18 # for permutation generation!

#SNP_DETAILS_GZ_FILENAME = "snp137.txt.gz" # Looks like this: 585     chr1    10258   10259   rs200940095     0       +       C       C       A/C     genomic single  unknown 0       0       near-gene-5     exact   1               1       GMI,    0 
# Tells you the actual letters. No header, so you have to guess what each column means. What a pain.

weirdHARsWrongAnnotSet = sets.Set() # List of HARs where we see something like: "snp is A/C" but then the actual sequence is a G or something. That is weird!


procCheck = subprocess.Popen(["egrep", "[	][-][	]", SNP_LETTERS_FILENAME], stdout=subprocess.PIPE) # specifically look for "tab hyphen tab" - which indicates a NEGATIVE STRAND thing still found in the file. We CANNOT OPERATE ON NEGATIVE STRANDS! You have to preprocess the file to switch the SNP locations to positive.
(out, err) = procCheck.communicate()


NUM_CHARS_THAT_MEANS_MATCH_PROBABLY_FOUND_SOMETHING = 10 # just sort of a ballpark "very small" number
if (len(out) > NUM_CHARS_THAT_MEANS_MATCH_PROBABLY_FOUND_SOMETHING):
    print "SNPs on negative strand: ", out
    print "Dang, super bad news!\n"
    print "The required input SNP file <" + SNP_LETTERS_FILENAME + "> had some NEGATIVE STRAND snps in it!"
    print "Unfortunately we don't know how to deal with those in this file. You need to"
    print "change all the SNPs so they are positive strand and flip them around accordingly."
    print "Example:   negative strand A/C would become positive strand T/G. Order does not matter."
    print " ---------- EXITING NOW ---------- "
    print " ---------- FIX THAT FILE! ---------- "
    sys.exit(1) # EXIT, FATAL ERROR!!!!!!!!!!!!!!!!!!!
    pass

if not os.path.exists(HAR_DIR):
    print "UH OH, the 'cut up hars' directory (OUT1_CUTUP_HARS) doesn't exist!! Run Hane's cut-up-HARs-to-remove-illegal-sequences script FIRST!"
    sys.exit(1)
    pass


def load_sequences(fasta_file):
    file1 = open(fasta_file, 'r')
    data = file1.readlines()
    
    str_data = ''.join(data)
    seq_data = filter(None, str_data.split('>'))

    hum_data = filter(None, seq_data[0].split('\n'))
    hum_seq = ''.join([x.strip() for x in hum_data[1:]])

    chimp_data = filter(None, seq_data[1].split('\n'))
    chimp_seq = ''.join([x.strip() for x in chimp_data[1:]])

    return (hum_seq, chimp_seq)


def strChange(string, loc, newChar):
    x = list(string)
    x[loc] = newChar
    return "".join(x)

def recursivelyMutate(theName, seq, changeHistory, remainSnps, remainLocs): # this is a RECURSIVE (!!) function
    if (len(remainSnps) == 0):
        return [ (seq, changeHistory) ] # nothing left to change, so just return the final sequence

    remainSnpsCOPY = copy.deepcopy(remainSnps) # copy this or else it totally destroys everything
    remainLocsCOPY = copy.deepcopy(remainLocs)
    snpList = remainSnpsCOPY.pop() # grab one of the snp options off the list, and we'll change the seq to have these snps
    snpLoc  = remainLocsCOPY.pop() # gotta know where the location in the HAR is
    savedResultsList = []
    for s in snpList: # once per snp option. e.g. this will happen 2x if the input is ["A","T"], or 3 times for ["A","C","G"]
        newChangeHistory = copy.deepcopy(changeHistory)
        newseq = str(seq) # copy the sequence
        
        if (snpLoc >= len(seq)):
            print "Looks like we have to SKIP the snp at location " + str(snpLoc) + " because it was GREATER THAN the bounds of the end of the sequences. That means it is in the NEXT split-up segment."
        elif snpLoc < 0:
            print "Looks like we have to SKIP the snp at location " + str(snpLoc) + " because it was LESS THAN THAN the bounds of this sequence. That means it was in the PREVIOUS split-up segment."
        else:
            # This SNP actually fell within the bounds of this split-up segment, so we should modify it.
	    # Let's make SUPER SURE that the base already at seq[snpLoc] is one of the ones that is in the annotation as a SNP
	    # So for example, if the first base is supposed to be either "A/G" according to the annotation, then seq[snpLoc] had better be
	    # either A or G. If it is a "T", then that probably means that something is seriously messed up with the annotation!
	    if (not (seq[snpLoc] in snpList)):
		if ((theName == "2xHAR.346") and (snpLoc == 52)):
		    print "This SNP is ok (manually verified). It is T in Chimp, and A/C in human."
		elif ((theName == "2xHAR.40") and (snpLoc == 4)):
		    print "This SNP is ok (manually verified). It is A in Chimp, and C/G in human."
		elif ((theName == "2xHAR.476") and (snpLoc == 23)):
		    print "This SNP is ok (manually verified). It is C in Chimp, and G/T in human."
		elif ((theName == "2xHAR.476") or (theName == "2xHAR.324") or (theName == "2xHAR.327") \
or (theName == "2xHAR.79") or (theName == "2xHAR.386") or (theName == "HAR64") \
or (theName == "2xHAR.102") or (theName == "2xHAR.62") or (theName == "2xHAR.214") or (theName == "2xHAR.451")):
		    print "This SNP is maybe ok, although I didn't check it!!!"
		else:
		    print "Whoa, possibly huge problem in sequence ", seq, "in", theName, "---the snp at location", snpLoc, "was supposed to be one of:", snpList, "but it was actually", seq[snpLoc]
		    print "There is probably some problem in your input annotation not matching up exactly! This may be a HUGE PROBLEM AND YOU NEED TO FIX YOUR ANNOTATION!!! You probably have an off-by-one error."
		    print "However, it is possible that we are actually looking at the CHIMP sequence and finding that for some (unlikely) reason, the human SNPs do not include the chimp sequence."
		    print "This should almost never happen, but it may happen a few times."
		    print "Note that this same function handles chimp AND human sequences with no differentiation."
		    global weirdHARsWrongAnnotSet
		    weirdHARsWrongAnnotSet.add(theName)
		    sys.exit(1) # <-- EXPLODE and quit!!!! This is super critical, or else you won't be able to find any problems
		    pass
	    
            if (seq[snpLoc] != s): # so we actually changed something!
                newChangeHistory.append((str(snpLoc+1)+""+ seq[snpLoc] +">"+s)) # note: snpLoc + 1 because it's indexed from 0 in python, but one in normal human being space!
                newseq = strChange(seq, snpLoc, s)
                pass

        m2 = recursivelyMutate(theName, newseq, newChangeHistory, remainSnpsCOPY, remainLocsCOPY)
        savedResultsList.extend(m2)
        pass
    return(savedResultsList)


def get_human_chimp_differences(hum_seq, chimp_seq):
    # Returns a TWO ELEMENT tuple, where the first is the difference array, second is the location array
    diff_array = [] # save the human/chimp differences, like so:  [ [A,C], [A,G], [T,C] }
    loc_array  = [] # save the LOCATION of the differences too!
    # first one is human!!!!!!
    if len(hum_seq) == len(chimp_seq):
        #diff_count = 0
        for i in range(0, len(hum_seq)):
            chimp_char = chimp_seq[i]
            hum_char   = hum_seq[i]
            if chimp_char != hum_char:
                #diff_count += 1
                diff_array.append( [ hum_char, chimp_char ] ) # first one is human!!!!
                loc_array.append( i )
                #print "DEBUG:", hum_char, chimp_char, i
                #print(','.join(map(str, [hum_char, chimp_char, i])) + '\n')
                pass
            pass
        pass
    else:
        print "UH OH LENGTH WAS NOT THE SAME! THIS IS A PROBLEM!!!!!!!!!!!!!!!!"
        sys.exit(1)
        pass

    return (diff_array, loc_array)


def get_seq_difference_count(seq1, seq2):
    diff_count = 0
    if len(seq1) == len(seq2):
        for i in range(0, len(seq1)):
            if seq2[i] != seq1[i]:
                diff_count += 1
                pass
            pass
        pass
    else:
        print "UH OH LENGTH WAS NOT THE SAME! THIS IS A PROBLEM!!!!!!!!!!!!!!!!"
        sys.exit(1)
        pass
    return diff_count


def stuff_to_fasta_format(har_name_stuff, baselinehuman, baselinechimp, snpList):
    # har_name_stuff is a string that should go at each fasta line
    # expects snpList to look something like this:
    # [('AAACCACTCAT', []), ('AAACCACTCATTTGTT', ['13 T -> G']), 
    #     ^             ^                                  ^
    #   the seq       changes from reference (a list)      note: this one only has ONE change (T->G at pos 13) but there can be tons of these
    # empty list = no changes from the reference
    fa = ""
    for item in snpList:
        # note: item is a TWO ELEMENT TUPLE!!!!!!!!!!!
        # first element: the sequence
        # second element: the ARRAY (list) of changes from the reference (possibly empty)
        # see above for examples
        seq = item[0]
        changeArray = item[1]
        numChanges = len(changeArray)

        numDiffsFromHumanReference = get_seq_difference_count(seq, baselinehuman)
        numDiffsFromChimpReference = get_seq_difference_count(seq, baselinechimp)


        ng  = len(re.findall("GGGG", seq))
        nr1 = len(re.findall("TCTAGA", seq))
        nr2 = len(re.findall("AGATCT", seq))

        if ((ng + nr1 + nr2) > 0):
            print "Seq is:" + seq
            print "How many GGGG:" + str(ng)
            print "How many TCTAGA:" + str(nr1)
            print "How many AGATCT:" + str(nr2)
            pass

        fa += ">" + har_name_stuff + "---" + str(numChanges) + "_modification(s)" + "---" + str(numDiffsFromHumanReference) + "_diff_human_ref" + "---" + str(numDiffsFromChimpReference) + "_diff_chimp_ref" + "---" + ":".join(changeArray) + "\n"
        fa += seq + "\n"
        pass
    return fa

# =====================================================

notFoundRSIDArray   = []
tooManyChangesArray = [] # save the har names that have too many things to change

rsidDict = {} # this will look like this:  rs10341  :  ["A", "C"]. Each rsID is associated with an array of possible snp values. Sometimes there can be MORE THAN TWO!!!
snpLettersFile = open(SNP_LETTERS_FILENAME, 'r')
snpData = snpLettersFile.readlines()
for line in snpData:
    if (len(line) <= 2): # string length... there is always an extra blank line at the end due to unix fun
        print "INFO: Ignoring an almost-blank line" # <-- ignore that line
        continue
    #print "OK" + line
    splitLine   = line.strip("\n").split("\t")
    theID       = splitLine[0]
    theLocation = int(splitLine[4]) # just for double checking
    theLetters  = splitLine[6]
    theLettersArray = theLetters.split("/")
    #print "Ok the letters array was" , theLettersArray
    rsidDict[ theID ] = {"loc":theLocation, "lettersArray":theLettersArray} # <-- make a tuple!!!
    #print rsidDict[theID]["loc"] # location
    #print rsidDict[theID]["lettersArray"] # snp 
    
    pass



#all_har_names = ["ABC123","SNAKES"]
all_har_paths = glob.glob(HAR_DIR + "/" + "*.fa") # os.listdir(HAR_DIR)

for har_path in all_har_paths:
    har_filename = os.path.basename(har_path)
    har_name = re.sub(".*_", "", re.sub(".fa", "", har_filename)) # turn "'new_chr17_2xHAR.123.fa'" into "2xHAR.123"
    (har_seq, chimp_seq)     = load_sequences(har_path)
    
    m = re.search("range_inclusive=(?P<STARTLOCATION>\d+)_", har_filename) # Match
    if (len(m.groups()) == 0):
        startCoordinateForThisHarStartingFromZero = 0
    else:
        startCoordinateForThisHarStartingFromZero = int(m.group("STARTLOCATION")) - 1 # subtract 1! go from 1-based indexing to 0-based. The startCoordinateForThisHarStartingFromZero tells us if this is a SPLIT UP har that didn't actually start from the HAR location in the annotation. There are certain "illegal" sequences in some HARs that require us to split up the HAR into two sub-hars. The SECOND (and onward) one of those won't match the annotation anymore, because it was split in half, so the start location is now incorrect. That's why we have to look for the literal text "range_inclusive=" some number, and then offset the HAR by the amount of the start of that range (minus one).
        pass
    
    print "The offset of this HAR was " + str(startCoordinateForThisHarStartingFromZero)

    print "----------------------------"
    print "HAR NAME: " + har_name
    
    # \\bhar2\\b <- the "\b" means "has to have something besides regular text here! (\b = word boundary) so now it won't find har231
    egrepHarStr = "\\b" + har_name + "\\b" # problem: searching for "har23" also finds "har230" and "har231" -- that's why we need \b, which requires a WORD BOUNDARY!
    proc    = subprocess.Popen(["egrep", egrepHarStr, HAR_COORDS_FILENAME], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()

    outSplit = out.strip("\n").strip("\r").split("\t")
    #print "program output:", out
    (egrepTerminalSTDOUT, egrepTerminalSTDERR) = subprocess.Popen(["egrep", egrepHarStr, SNP_PER_HAR_FILENAME], stdout=subprocess.PIPE).communicate()
    snpChrs   = []
    snpStarts = []
    snpEnds   = []
    snpNames  = [] # <-- these are rsID identifiers. You can get more info about each snp from the 'rsidDict'
    for line in egrepTerminalSTDOUT.split("\n"):
        # Getting all the SNPs out of this file...
        if (len(line) <= 2): # string length... there is always an extra blank line at the end due to unix fun
            print "INFO: Ignoring an almost-blank line" # <-- ignore that line
            continue
        x = line.strip("\r").split("\t")

        theRSID = x[3]

        if (theRSID in snpNames):
            # dangit, this is a duplicate rsID in the annotation! Skip it
            print "Skipping the DUPLICATE rsID in the annotation <" + theRSID + "> (we already added it)"
            continue

        snpChrs.append(x[0]) # chrom name
        snpStarts.append(int(x[1])) # number
        snpEnds.append(int(x[2])) # number
        snpNames.append(theRSID)
        if ((not x[5] in har_name) and (not har_name in x[5])): # har name, for sanity checking
            print "Here was the name in x5: ", x[5]
            print "Here was the har name:", har_name
            print "HAR NAME DID NOT MATCH! THIS IS A SERIOUS ERROR"
            sys.exit(1)
            pass
        pass

    #print "Here are the snps associated with this HAR: ", snpNames
    #sys.exit(1)

    snp_rsid_array         = []
    snp_diff_array         = []
    snp_loc_array_absolute = []

    for rs in snpNames:
        #print rs
        #print rsidDict[rs]
        if (rs not in rsidDict):
            print "HEY!!!! The rsID was not found in the dictionary!!! That means that this rsID was IN your original data file, but it was NOT in the 'which letter is this SNP' file, so we ignored it!"
            notFoundRSIDArray.append(rs)
            continue # welp, forget about this one, go to the next one
        theDifferenceArray = rsidDict[rs]["lettersArray"]
        theDifferenceLoc   = rsidDict[rs]["loc"]
        snp_rsid_array.append(rs) # save the rsID
        snp_diff_array.append(theDifferenceArray)
        snp_loc_array_absolute.append(theDifferenceLoc)
        pass

    print "Outsplit was: " , outSplit, " (if this is empty, it means we didn't find the HAR in the annotation file, which causes a CRASH)"
    if (len(outSplit) == 1):
        print "SERIOUS ERROR::::::::::::::::::::::::::"
        print "outSplit (which is the output of the 'egrep' command"
        print "did not find any entries for this HAR. Maybe it was a testing har that is not real, or maybe your annotation is in the wrong place?"
        print "A CRASH IS GOING TO HAPPEN IN A FEW LINES HERE"
        pass

    har_chr   = outSplit[0]
    har_real_start = int(outSplit[1]) # <-- the "real" start and end are the locations for the actual har, NOT our split-up hars
    har_real_end   = int(outSplit[2])

    har_substart = har_real_start + startCoordinateForThisHarStartingFromZero # This accounts for the possibility of SPLITTING the har up to avoid the illegal sequences
    har_subend   = (har_substart + len(har_seq))

    print "From the har coordinates file (NOT the filename!) <" + har_filename + "> we extracted har <" + har_name + "> at real non-split-up position " + har_chr + ":" + str(har_real_start) + "-" + str(har_real_end)
    print "From the filename, however, we found that the split-up har was actually at " + har_chr + ":" + str(har_substart) + "-" + str(har_subend)
    print "Human seq was: " + har_seq
    print "Chimp seq was: " + chimp_seq
    # get a list of the SNPs that are associated with this HAR

    (chimp_human_diff_array, chimp_human_loc_array) = get_human_chimp_differences(har_seq, chimp_seq)

    # chimp_human_loc_array is the RELATIVE positions in the har where the human-chimp differences are. Example: 0 would mean "the first base is different in human/chimp".
    #   * It is similar to the variable "chimp_human_diff_loc_relative"

    num_diff_human_chimp = len(chimp_human_loc_array) # get_seq_difference_count(har_seq, chimp_seq)
    print "This many differences: " + str(num_diff_human_chimp)
    
    #snps_this_har                            = [["A","B"], ["C","D"],["E","F","G","H"]] # get them from the SNP_PER_HAR_FILENAME
    chimp_human_diff_loc_relative = map(lambda x: x - har_substart, chimp_human_loc_array) # <-- this has to be a SMALL NUMBER!!!!!!!!!!!!! Like, 0 = start of HAR
    snp_loc_array_relative        = map(lambda x: x - har_substart, snp_loc_array_absolute) # <-- this has to be a SMALL NUMBER!!!!!!!!!!!!! Like, 0 = start of HAR.

    #print chimp_human_diff_array
    #print chimp_human_loc_array
    #print snp_diff_array
    #print snp_loc_array_absolute

    #if (snp_loc_array_absolute != snpStarts):
    #    print "THIS IS A SERIOUS PROBLEM -- one file gave the snpStarts as one thing, but another gave them as another thing! Debug this, because it is a huge possible problem!!!!!!!!!!!!!!!"
    #    print "Note that we are ignoring it in the cases when the snps have MULTIPLE locations listed"
    #    print snp_loc_array_absolute
    #    print snpStarts
    #    sys.exit(1)
    #    pass

    print "SNP CHANGE LIST on", har_chr, ":", snp_diff_array
    print "ABSOLUTE SNP LOCATION (hg19) on", har_chr, ":", snp_loc_array_absolute
    print "RELATIVE SNP LOCATION:", snp_loc_array_relative

    print "CHIMP/HUMAN DIFF:", chimp_human_diff_array
    print "CHIMP/HUMAN LOCATIONS (relative to seq. start):", chimp_human_loc_array

    #sys.exit(1)
    
    #snpifiedList   = recursivelyMutate(har_seq, [], snp_diff_array, snp_loc_array_relative)
    #chimpifiedList = recursivelyMutate(har_seq, [], chimp_human_diff_array, chimp_human_loc_array)
    #print (snp_diff_array + chimp_human_diff_array)
    #print (snp_loc_array_relative + chimp_human_loc_array)

    both_diff_array = (snp_diff_array + chimp_human_diff_array)
    both_loc_array  = (snp_loc_array_relative + chimp_human_loc_array)

    unique_diff_array = []
    unique_loc_array  = []

    # ============================ This is all a giant function to only extract the UNIQUE changes ================
    seenSeqAndLoc = {}
    for i in range(0, len(both_diff_array)):
        thisDiff = both_diff_array[i]
        thisLoc  = both_loc_array[i]

        if ('*' in thisDiff):
            print "We are totally skipping this INDEL!"
            continue

        if (thisLoc not in seenSeqAndLoc):
            # we've never seen any mutation at THIS location, so add this one!
            unique_loc_array.append(thisLoc)
            unique_diff_array.append(thisDiff)
            seenSeqAndLoc[thisLoc] = [] # new empty array! We're going to save all the mutations we see at this point here. So if we see a [A,C], we will APPEND it to this list. If we later see an [A,T], we'll append that too!
            seenSeqAndLoc[thisLoc].append(thisDiff)
        else:
            # we have seen something at this location before, let's see if
            # the DIFFERENCE is the same difference we saw before. i.e. is it [A,C] and [A,C] again, or is it something new like [A,C] and [A,T]
            haveWeSeenItBefore = False
            for x in seenSeqAndLoc[thisLoc]:
                if (set(x) == set(thisDiff)):
                    haveWeSeenItBefore = True
                    pass
                pass
            
            if (haveWeSeenItBefore):
                # yep, we've seen this exact one before! don't add this change
                print "Hey we saw a duplicate change at " + str(thisLoc) + ", so we are not adding it again"
            else:
                # looks like it was a different type of mutation, so add it!
                seenSeqAndLoc[thisLoc].append(thisDiff)
                unique_loc_array.append(thisLoc)
                unique_diff_array.append(thisDiff)
                pass
            pass
        pass


    print "There were this many human SNP differences: " + str(len(snp_diff_array))
    print "There were this many human/chimp differences: " + str(len(chimp_human_diff_array))

    print "There were this many OVERALL differences: " + str(len(unique_diff_array))

    if (len(unique_diff_array) > MAX_ALLOWED_DIFFERENCES):
        rsjoin = ",".join(map(str, snp_rsid_array))
        hjoin = ",".join(map(str, snp_diff_array))
        hloc  = ",".join(map(str, snp_loc_array_relative))

        cjoin = ",".join(map(str, chimp_human_diff_array))
        cloc  = ",".join(map(str, chimp_human_loc_array))

        sss = har_name \
            + " had " \
            + str(len(snp_diff_array)) \
            + " snp changes" \
            + " (" + rsjoin + ") " \
            + " (" + hjoin + ") " \
            + " (" + hloc + ") " \
            + "and " \
            + str(len(chimp_human_diff_array)) \
            + " human/chimp changes" \
            + " (" + cjoin + ") " \
            + " (" + cloc + ") "

        tooManyChangesArray.append(sss)
        continue # skip it entirely, too many changes!!!

    bothifiedList   = recursivelyMutate(har_name, har_seq, [], unique_diff_array, unique_loc_array) # note: this is both human AND chimp sequences
    
    #print chimp_human_diff_loc_relative

    print "The permutations of all SNP differences are", bothifiedList
    print "The permutations of all chimp-human differences (", str(len(bothifiedList)), "in all) are", bothifiedList

    har_stuff = har_filename + "---" + har_name + "---" + str(har_chr) + ":" + str(har_real_start) + "-" + str(har_real_end) + "---" + "split_location_at_" + str(har_chr) + ":" + str(har_substart) + "-" + str(har_subend)
    # =====================================================================================

    # ===================== kick out any duplicates that have inexplicably STILL remained here!
    uniqueList = []
    seenSeq     = {}    
    for item in bothifiedList:
       seq = item[0]
       if (seq not in seenSeq):
           uniqueList.append(item)
       else:
           print "We already saw '" + seq + "', so we aren't adding it again!"
           pass
       seenSeq[seq] = True
       pass
    print "Started with " + str(len(bothifiedList)) + " items, and ended up with " + str(len(uniqueList)) + " after removing duplicates."
    # ============================================================

    outFastaString = stuff_to_fasta_format(har_stuff, har_seq, chimp_seq, uniqueList)
    
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
        pass

    outFilename = "processed." + har_filename
    outPath     = OUT_DIR + "/" + outFilename
    
    writeFile = open(outPath, 'w')
    writeFile.write(outFastaString)
    writeFile.close()

    print "Here are the RSIDs that were listed in one file but not found in the 'which letter' file:"
    print "\n".join(notFoundRSIDArray)

    print "Here are the HARs that had too many SNPs or perhaps human/chimp changes, where 'too many' is MORE THAN " + str(MAX_ALLOWED_DIFFERENCES) + ":"
    print "\n".join(tooManyChangesArray)

    print "Here are the " + str(len(weirdHARsWrongAnnotSet)) + " HARs with 'weird' SNPs where the annotation said something like 'snp is A/G' but the actual genome was a T or C (e.g. not A/G): "
    print ", ".join(weirdHARsWrongAnnotSet) + "\n"

    pass


# 38049 permutations for all MODIFIED SPLIT UP HARS in : /home/hane/ENHANCERFINDER_HARs/MODIFIED_all_permutations.fa

#################################################################################################################
# SIXTH STEP. COVER ALL PERMUTATIONS FOR ENHANCERFINDER HARS
#################################################################################################################

#to grab all permutations for EnhancerFinder HARs, see 

cut -f 4 239_EnhancerFinder_HARs.bed | perl -p -e 's/;/\n/' | perl -n -e "chomp($_); print '|_'; print $_; print '_'" | sed 's/|//' > grepcommand.txt

# MAKES THE HORRIBLE GREP COMMAND

x=`cat grepcommand.txt`
x="($x)"

# add the parens and make x into a variable in bash/sh

now you can type echo $x to see x

egrep -A 1 "$x" MODIFIED_all_permutations.fa | grep -v '^--' > 239_EnhancerFinder_HAR_permutations.fa

# counts for each guy, assumes '_' split
cat interesting.fa | grep '^>' | cut -d '_' -f 3 | less -S | sort | uniq -c > 239_EnhancerFinder_HAR_permutations.txt

#################################################################################################################

# make a list of EnhancerFinder HARs that all permutations are covered for
import itertools

directory = '/home/hane/ENHANCERFINDER_HARs/'

file1 = open(directory + '239_EnhancerFinder_HARs.bed', 'rb')
file2 = open(directory + 'EnhancerFinder_list.txt', 'wb')

EnhancerFinder_aliases = []
EnhancerFinder_HARs = []

for line in file1:
	data = line.strip('\r\n').split('\t')
	chrom = data[0]
	start = data[1]
	end = data[2]
	name =data[3]
	if ';' in name:
		aliases = name.split(';')
		EnhancerFinder_aliases.append(aliases)
	else:
		EnhancerFinder_HARs.append(name)

merged_EnhancerFinders_HARs = list(itertools.chain.from_iterable(EnhancerFinder_aliases))

HAR_aliases = []
for HAR in merged_EnhancerFinders_HARs:
	if 'HAR' in HAR:
		HAR_aliases.append(HAR)
merged_lists = HAR_aliases + EnhancerFinder_HARs

for x in merged_lists:
	file2.write(x + '\n')
file2.close() 

#################################################################################################################

# make a BED file of HARs that are not EnhancerFinder predictions for which all permutations will be covered 

import collections

hars_directory = '/Users/haneryu/Documents/Pollard_Lab_2013/Prioritizing_permutations/'
directory = '/Users/haneryu/Documents/Pollard_Lab_2013/'

file1 = open(directory + 'EnhancerFinder_list.txt', 'rb')
file2 = open(hars_directory + 'sorted_hars_hg19.bed', 'rb')
file3 = open(directory + 'MODIFIED_non_EnhancerFinder_HARs.txt', 'wb')
file4 = open(directory + 'MODIFIED_non_EnhancerFinder_HARs.bed' ,'wb')

EnhancerFinder_list = sorted(set(list([])))
for line in file1:
	HAR_name = line.strip('\n')
	EnhancerFinder_list.append(HAR_name)
print 'EnhancerFinder_list :' , len(EnhancerFinder_list)

all_HAR_dictionary = {}
all_HARs = []
all_HAR_lines = []
for line in file2:
	all_HAR_lines.append(line)
	data = line.strip('\r\n').split('\t')
	chrom = data[0]
	start = data[1]
	end = data[2]
	HAR_identifier = data[3]
	if ';' in HAR_identifier:
		HAR_aliases = HAR_identifier.split(';')
		HAR_identifier = HAR_aliases[0]
		#print first_alias
		#second_alias = HAR_aliases[1]
		#print second_alias
		all_HARs.append(HAR_identifier)
		#all_HARs.append(second_alias)
	else:
		#print HAR_identifier
		all_HARs.append(HAR_identifier)
	all_HAR_dictionary[HAR_identifier] = line
print 'all_HARs:', len(all_HARs)
print 'all_HAR_dictionary:' , len(all_HAR_dictionary)

EnhancerFinder_list_multiset = collections.Counter(EnhancerFinder_list)
all_HARs_multiset = collections.Counter(all_HARs)

non_EnhancerFinder_HARs = sorted(set(list((all_HARs_multiset - EnhancerFinder_list_multiset).elements())))

print 'non_EnhancerFinder_HARs:' , len(non_EnhancerFinder_HARs)

for HAR in non_EnhancerFinder_HARs:
	BEDfile_line = all_HAR_dictionary[HAR]
	print BEDfile_line
	file4.write(BEDfile_line)
file4.close()

#################################################################################################################
# SEVENTH STEP. PERFORM MATCH CALL AS A STEP IN PRIORITIZING PERMUTATIONS BASED ON TFBS MISMATCHES
#################################################################################################################

# MAKE A nonEnhancerFinder_HARs.fa file using the comments.txt commands 


# append HAR name in all identifiers in split up HARs:

# add HAR name to identifier in split up HAR FASTA files

import os

directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/'
output_directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/2xHARs/'

def load_sequences(fasta_file):
		file1 = open(directory + fasta_file, 'rb')
		data = file1.readlines()
		str_data = ''.join(data)
		seq_data = filter(None, str_data.split('>'))

		hum_data = filter(None, seq_data[0].split('\n'))
		hum_identifier = hum_data[0]
		hum_seq = ''.join([x.strip() for x in hum_data[1:]])

		chimp_data = filter(None, seq_data[1].split('\n'))
		chimp_identifier = chimp_data[0]
		chimp_seq = ''.join([x.strip() for x in chimp_data[1:]])

		return hum_seq, hum_identifier, chimp_seq, chimp_identifier

def grab_HAR_name(fasta_file):
	if fasta_file.endswith('.fa'):
		fasta_file_name = fasta_file.split('_')
		HAR_name = fasta_file_name[2]
		return HAR_name, fasta_file


fasta_files = os.listdir(directory)

for fasta_file in fasta_files:
	if fasta_file.endswith('.fa'):
		if fasta_file.startswith('split_up_2x'):
			hum_seq, hum_identifier, chimp_seq, chimp_identifier = load_sequences(fasta_file)
			HAR_name, fasta_file = grab_HAR_name(fasta_file)
			print HAR_name
			print fasta_file
			output_file_name = fasta_file
			output_file = open(output_directory + fasta_file, 'wb')
			#print hum_identifier
			#print hum_seq
			#print chimp_identifier
			#print chimp_seq
			parsed_hum_identifier = hum_identifier.split('__')
			parsed_chimp_identifier = chimp_identifier.split('__')
			hum_species = parsed_hum_identifier[0]
			chimp_species = parsed_chimp_identifier[0]
			hum_cut = parsed_hum_identifier[1]
			chimp_cut = parsed_chimp_identifier[1]
			hum_length = parsed_hum_identifier[2]
			chimp_length = parsed_chimp_identifier[2]
			hum_range = parsed_hum_identifier[3]
			chimp_range = parsed_chimp_identifier[3]
			#print hum_identifier
			#print hum_species, HAR_name, hum_cut, hum_length, hum_range
			#print chimp_identifier
			#print chimp_species, HAR_name, chimp_cut, chimp_length, chimp_range
			new_hum_identifier = str( '>' + hum_species + '__' + HAR_name + '__' + hum_cut + '__' + hum_length + '__' + hum_range)
			new_chimp_identifier = str( '>' + chimp_species + '__' + HAR_name + '__' + chimp_cut + '__' + chimp_length + '__' + chimp_range)
			print new_hum_identifier
			print new_chimp_identifier
			output_file.write(new_hum_identifier + '\n')
			output_file.write(hum_seq + '\n')
			output_file.write(new_chimp_identifier + '\n')
			output_file.write(chimp_seq + '\n')
output_file.close()

#################################################################################################################

#GRAB ALL HG19 and PT2 SEQUENCES PER SPLIT UP HAR from /home/hane/MODIFIED_SPLIT_UP_HARs/MODIFIED_SPLIT_UP_HARs_1067.fa:

# all in this directory in lighthouse: /home/hane/nonEnhancerFinder_HARs

cut -f 4 MODIFIED_non_EnhancerFinder_HARs.bed | perl -p -e 's/;/\n/' | perl -n -e "chomp($_); print '|_'; print $_; print '_'" | sed 's/|//' > grepcommand.txt
# MAKES THE HORRIBLE GREP COMMAND

x=`cat grepcommand.txt`
x="($x)"

# add the parens and make x into a variable in bash/sh

now you can type echo $x to see x

egrep -A 1 "$x" MODIFIED_SPLIT_UP_HARs_1067.fa | grep -v '^--' > MODIFIED_non_EnhancerFinder_HARs.fa

# counts for each guy, assumes '_' split
cat MODIFIED_non_EnhancerFinder_HARs.fa | grep '^>' | cut -d '_' -f 3 | less -S | sort | uniq -c > MODIFIED_non_EnhancerFinder_HARs.txt

#file is here : /home/hane/nonEnhancerFinder_HARs/MODIFIED_non_EnhancerFinder_HARs.fa

#################################################################################################################

# match call command: 
# /data/transfac/1_transfac/transfac_2012.4/match/bin/match /data/transfac/1_transfac/transfac_2012.4/match/data/matrix.dat /home/hane/nonEnhancerFinder_HARs/MODIFIED_non_EnhancerFinder_HARs.fa /home/hane/MODIFIED_TFBS_MISMATCH_FILES/nonEnhancerFinder_TFBS_output.match /data/transfac/1_transfac/transfac_2012.4/match/data/prfs/vertebrate_non_redundant_minSUM.prf

# ALL TFBS OUTPUT FILES GO TO: /home/hane/MODIFIED_TFBS_MISMATCH_FILES

#################################################################################################################
# SEVENTH STEP. IDENTIFY ALL SPLIT UP HARS THAT CONTAIN TFBS MISMATCHES
#################################################################################################################

# need to create one-based nonEnhancerFinder BED file
# parse split up HAR fasta files to generate a BED file that has a fake chrom, start site = 0, end site = length -1 , identifier ---> so Tara can run her TFBS-HAR mutation locater program

import os

directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/'
output_BEDfile = open(directory + 'one_based_nonEnhancerFinder_HARs.bed', 'wb')

split_up_HARs = os.listdir(directory)

for split_up_HAR in split_up_HARs:
	if split_up_HAR.endswith('.fa'):
		file1 = open(directory + '/' + split_up_HAR, 'rb')
		data = file1.readlines()
		str_data = ''.join(data)
		#print str_data
		seq_data = filter(None, str_data.split('>'))
		#print seq_data
		hum_data = filter(None, seq_data[0].split('\n'))
		hum_identifier = hum_data[0]
		split_HAR_info = hum_identifier.split('_')
		hum_cut_number = split_HAR_info[4]
		hum_length = split_HAR_info[8]
		split_HAR_length = hum_length.split('=')
		oligo_length = split_HAR_length[1]
		end_site = int(oligo_length) + 1

		#print hum_data
		hum_seq = ''.join([x.strip() for x in hum_data[1:]])
		#print hum_seq
		chimp_data = filter(None, seq_data[1].split('\n'))
		#print chimp_data
		chimp_identifier = chimp_data[0]
		chimp_split_HAR_info = chimp_identifier.split('_')
		chimp_cut_number = chimp_split_HAR_info[2]
		chimp_length = chimp_split_HAR_info[6]

		chimp_seq = ''.join([x.strip() for x in chimp_data[1:]])
		#print chimp_seq
		print hum_identifier, hum_seq, hum_cut_number, hum_length, oligo_length, end_site
		print chimp_identifier, chimp_seq, chimp_cut_number, chimp_length
		hum_row_to_write = ['chr_unknown', '1', int(end_site), hum_identifier ]
		print hum_row_to_write
		output_BEDfile.write('\t'.join(map(str, hum_row_to_write)) + '\n')

		chimp_row_to_write = ['chr_unknown', '1', int(end_site), chimp_identifier ]
		print chimp_row_to_write

		output_BEDfile.write('\t'.join(map(str, chimp_row_to_write)) + '\n')
output_BEDfile.close()

#################################################################################################################

# need to create a modified non EnhancerFinder hg19-pt2 diffs BED file

#FIRST CREATE A DIFFERENCES FILE DIRECTORY:

import os

def load_sequences(fasta_file):
		file1 = open(fasta_file, 'rb')
		data = file1.readlines()
		str_data = ''.join(data)
		seq_data = filter(None, str_data.split('>'))

		hum_data = filter(None, seq_data[0].split('\n'))
		hum_seq = ''.join([x.strip() for x in hum_data[1:]])

		chimp_data = filter(None, seq_data[1].split('\n'))
		chimp_seq = ''.join([x.strip() for x in chimp_data[1:]])

		return hum_seq, chimp_seq


def find_human_chimp_differences(hum_seq, chimp_seq, out_file):
		if len(hum_seq) == len(chimp_seq):
			diff_count = 0
		for i in range(0, len(hum_seq)):
			chimp_char = chimp_seq[i]
			hum_char = hum_seq[i]
			if chimp_char != hum_char:
				diff_count += 1
				print hum_char, chimp_char, i
				out_file.write('\t'.join(map(str, [hum_char, chimp_char, i])) + '\n')

directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/'
diff_directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/DIFFS/'

fasta_files = os.listdir(directory)
for fasta_file in fasta_files:
		#print fasta_file
		if fasta_file.endswith('.fa'):
			outfile_name = fasta_file.replace('.fa', '.diff')
			out_file = open(diff_directory + outfile_name, 'wb')
			hum_seq, chimp_seq = load_sequences(directory + '/' + fasta_file)
			#find_human_chimp_differences(hum_seq, chimp_seq, out_file)               
			find_human_chimp_differences(hum_seq, chimp_seq, out_file)
			#out_file.write('\t'.join(map(str, [hum_char, chimp_char, i])) + '\n')
		out_file.close()       

# total number of fixed human-chimp differences = 2817

#################################################################################################################

# take hg19_pt2 differences per split up HAR and put them in one large BED file with the locations relative to the start site of the split up HAR

import os

directory = '/home/hane/MODIFIED_SPLIT_UP_HARs/DIFFS/'
outfile = open(directory + 'non_EnhancerFinder_HARs_hg19_pt2_diffs.bed', 'wb')

def load_identifier(human_chimp_difference):
	diff_name = human_chimp_difference.strip('.diff').split('_')
	#print diff_name
	HAR_name = diff_name[2]
	space = diff_name[3]
	#print HAR_name
	cut_number = diff_name[4]
	of = diff_name[5]
	cut_total = diff_name[6]
	second_space = diff_name[7]
	#print cut_number
	length = diff_name[8]
	#print length
	range_inclusive = diff_name[9]
	inclusion_number = diff_name[10]
	to = diff_name[11]
	inclusion_end = diff_name[12]
	end = diff_name[13]
	new_name = [HAR_name, space, cut_number, of, cut_total, second_space, length, range_inclusive, inclusion_number, to, inclusion_end, end]
	new_identifier = '_'.join(map(str, new_name))
	return new_identifier, HAR_name, cut_number, length


human_chimp_differences = os.listdir(directory)

for human_chimp_difference in human_chimp_differences:
	if human_chimp_difference.endswith('.diff'):
		new_identifier, HAR_name, cut_number, length = load_identifier(human_chimp_difference)
		print new_identifier
		print HAR_name, cut_number, length
		file1 = open(directory + human_chimp_difference, 'rb')
		for line in file1:
			data = line.strip('\n').split('\t')
			human_char = data[0]
			chimp_char = data[1]
			diff_location = data[2]
			one_based_location = int(diff_location) + 1
			print human_char, chimp_char, diff_location, one_based_location
			hum_species = 'hg19__'
			hum_identifier = ''.join(map(str, [hum_species, new_identifier]))
			chimp_species = 'panTro2__'
			chimp_identifier = ''.join(map(str, [chimp_species, new_identifier]))
			hum_row_to_write = [one_based_location, hum_identifier]
			outfile.write('\t'.join(map(str, hum_row_to_write)) + '\n')
			print hum_row_to_write
			chimp_row_to_write = [one_based_location, chimp_identifier]
			print chimp_row_to_write
			outfile.write('\t'.join(map(str, chimp_row_to_write)) + '\n')
outfile.close()

#HG19-PT2 DIFF FILE IS HERE: /home/hane/MODIFIED_SPLIT_UP_HARs/DIFFS/non_EnhancerFinder_HARs_hg19_pt2_diffs.bed

#################################################################################################################
# EIGHTH STEP. IDENTIFY ALL SPLIT UP HARS THAT CONTAIN TFBS MISMATCHES
#################################################################################################################

# TFBS_output.match FILE IS HERE: /home/hane/MODIFIED_TFBS_MISMATCH_FILES/nonEnhancerFinder_TFBS_output.match 
# ONE-BASED HAR BED FILE IS HERE: /home/hane/MODIFIED_TFBS_MISMATCH_FILES/one_based_nonEnhancerFinder_HARs.bed
# HG19-PT2 DIFF FILE IS HERE: /home/hane/MODIFIED_TFBS_MISMATCH_FILES/non_EnhancerFinder_HARs_hg19_pt2_diffs.bed

# RUN TARA's PROGRAM HERE: /home/tara/co/scripts

python runTFModels.py /home/hane/MODIFIED_TFBS_MISMATCH_FILES/one_based_nonEnhancerFinder_HARs.bed /home/hane/MODIFIED_TFBS_MISMATCH_FILES/nonEnhancerFinder_TFBS_output.match  2 /home/hane/MODIFIED_TFBS_MISMATCH_FILES/non_EnhancerFinder_HARs_hg19_pt2_diffs.bed > /home/hane/MODIFIED_TFBS_MISMATCH_FILES/TFBS_output.txt

# OUTPUT FILE: /home/hane/MODIFIED_TFBS_MISMATCH_FILES/TFBS_output.txt

#################################################################################################################

#parse the TFBS output.txt file to get two lists of HARs in which TFBS motifs only occur in the human OR the chimp sequences given
import collections

directory = '/home/hane/MODIFIED_TFBS_MISMATCH_FILES'

file1 = open(directory + '/' + 'TFBS_output.txt', 'rb')
file2 = open(directory + '/' + 'TFBS_mismatch_hg19.txt', 'wb')
file3 = open(directory + '/' + 'TFBS_mismatch_pt2.txt', 'wb')

hg19_KEYS = []
pt2_KEYS = []
hg19_MOTIFS = [] 
pt2_MOTIFS = []
hg19_oligos = []
pt2_oligos = []

#my_dictionary = {}
hg19_dictionary = {}
pt2_dictionary = {}


for line in file1:
	data = line.strip('\n').split('\t')
	#print data
	mismatch_position = data[0]
	motif_start = data[1]
	motif_end = data[2]
	strand_direction = data[3]
	motif_length = data[4]
	motif_seq = data[5]
	motif_tfname = data[6]
	region_length = data[7]
 	region_name = data[8]
 	region_name_parts = region_name.split('__')
 	print region_name_parts
 	species = region_name_parts[0] 
 	print species
 	HAR = region_name_parts[1]
	print HAR
	cut_number = region_name_parts[2]
	one_based_oligo_length = region_name_parts[3]
	range_inclusive = region_name_parts[4]
	#KEY = str(mismatch_position + '__' +  HAR + '__' + cut_number + '__' + one_based_oligo_length + '__' + range_inclusive )
	KEY = str( HAR + '__' + cut_number + '__' + one_based_oligo_length + '__' + range_inclusive )
	#print KEY
	MOTIF_IDENTIFIER = str(motif_tfname + '__' + motif_start + '__' + motif_end + '__' + motif_length)
	#print MOTIF_IDENTIFIER
	if 'hg19' in species:
		#print data
		hg19_KEYS.append(KEY)
		hg19_MOTIFS.append(MOTIF_IDENTIFIER)
		hg19_oligos.append(region_name)
		hg19_dictionary[KEY] = region_name
		#hg19_dictionary[KEY] = MOTIF_IDENTIFIER
		pass
	if 'panTro2' in species:
		pt2_KEYS.append(KEY)
		pt2_MOTIFS.append(MOTIF_IDENTIFIER)
		pt2_oligos.append(region_name)
		#pt2_dictionary[KEY] = region_name
		pt2_dictionary[KEY] = MOTIF_IDENTIFIER
		pass
#print hg19_KEYS
print 'hg19_KEYS LIST HAS' , len(hg19_KEYS) , 'ELEMENTS'
#print hg19_MOTIFS
print 'hg19_MOTIFS LIST HAS' , len(hg19_MOTIFS) , 'ELEMENTS'
#print hg19_dictionary
#print pt2_KEYS
print 'pt2_KEYS LIST HAS' , len(pt2_KEYS) , 'ELEMENTS'
#print pt2_MOTIFS
print 'pt2_MOTIFS LIST HAS' , len(pt2_MOTIFS) , 'ELEMENTS'
#print pt2_dictionary

hg19_HAR_identifier_multiset = collections.Counter(hg19_KEYS)
pt2_HAR_identifier_multiset = collections.Counter(pt2_KEYS)

non_mismatch_TFBS_HAR_identifiers = list((hg19_HAR_identifier_multiset & pt2_HAR_identifier_multiset).elements())
print 'number of HARs that have a TFBS motifs that are unbroken in both human and chimp:', len(non_mismatch_TFBS_HAR_identifiers)
#print non_mismatch_TFBS_HAR_identifiers

hg19_HAR_identifier_remainder = sorted(set(list((hg19_HAR_identifier_multiset- pt2_HAR_identifier_multiset).elements())))
pt2_HAR_identifier_remainder = sorted(set(list((pt2_HAR_identifier_multiset - hg19_HAR_identifier_multiset).elements())))

print 'number of hg19 HAR sequences that have a TFBS motif not in pt2:', len(hg19_HAR_identifier_remainder) # just in human
#print hg19_HAR_identifier_remainder

print 'number of pt2 HAR sequences that have a TFBS motif not in hg19:', len(pt2_HAR_identifier_remainder) # just in chimp
#print pt2_HAR_identifier_remainder


for har in hg19_HAR_identifier_remainder:
	#pt2_dictionary[har] = MOTIF_IDENTIFIER
	#har, MOTIF_IDENTIFIER
	#print har, MOTIF_IDENTIFIER
	file2.write( har +'\n')
	#file3.write('\t'.join(map(str, [har, MOTIF_IDENTIFIER]))+'\n')
file2.close()

for har in pt2_HAR_identifier_remainder:
	#pt2_dictionary[har] = MOTIF_IDENTIFIER
	#har, MOTIF_IDENTIFIER
	#print har, MOTIF_IDENTIFIER
	file3.write( har +'\n')
	#file3.write('\t'.join(map(str, [har, MOTIF_IDENTIFIER]))+'\n')
file3.close()

# hg19_KEYS LIST HAS 14180 ELEMENTS
# hg19_MOTIFS LIST HAS 14180 ELEMENTS
# pt2_KEYS LIST HAS 14833 ELEMENTS
# pt2_MOTIFS LIST HAS 14833 ELEMENTS
# number of HARs that have a TFBS motifs that are unbroken in both human and chimp: 12,681
# number of hg19 HAR sequences that have a TFBS motif not in pt2: 240
# number of pt2 HAR sequences that have a TFBS motif not in hg19: 298

#################################################################################################################

# DID THE FOLLOWING TO EXTRACT PERMUTATION SEQUENCES FOR ONLY TFBS MISMATCH CONTAINING HARS in: /home/hane/MODIFIED_TFBS_MISMATCH_FILES

x=`cat TFBS_mismatch_hg19_AND_pt2.txt`
x="$x"
egrep -A 1 "$x" MODIFIED_all_permutations.fa | grep -v '^--' > TFBS_mismatch_nonEnhancerFinder_HARs.fa



# run through TFBS_mismatch_nonEnhancerFinder_HAR_permutations.fa file and generate output files that only have human and chimp and single modification sequences

import sys
import re


directory = '/home/hane/MODIFIED_TFBS_MISMATCH_FILES/'
file1 = open(directory + 'TFBS_mismatch_nonEnhancerFinder_HARs.fa', 'rb')

def ReadFasta(fh):
	''' Reads fasta file and returns alignment as a dict.
	
	'''

	ali = dict()
	header = ''
	seq = ''

	for line in fh:
		line = line.rstrip('\n\r')
		if line == '':
			continue
		if line.startswith('>'):
			if header != '':
				ali[header] =  re.sub('\s+', '',seq)
				seq = ''
			header = line
			continue
		seq += line
	ali[header] =  re.sub('\s+', '',seq)
	return ali

x = ReadFasta(file1) 
#print x


hum_refs = []
chimp_refs = []
single_modifications = []
double_modifications = []
triple_modifications = []

for identifer in x.keys():
	oligo_descriptor = str(identifer)
	parsed_identifer = identifer.split('--')
	#print parsed_identifer
	modification = parsed_identifer[4]
	#print modification
	hum_ref = parsed_identifer[5]
	chimp_ref = parsed_identifer[6]
	if '-0_modification(s)' in modification:
		seq = x[oligo_descriptor]
		hum_refs.append(oligo_descriptor)
		#print oligo_descriptor
		#print seq
		# print x
		# print 'hum_ref'
		# print modification
	if '-1_modification(s)' in modification:
		seq = x[oligo_descriptor]
		#print oligo_descriptor
		#print seq
		single_modifications.append(oligo_descriptor)
		pass
		# print x
		# print 'single_modification'
		# print modification
	if '-0_diff_chimp_ref' in chimp_ref:
		seq = x[oligo_descriptor]
		#print oligo_descriptor
		#print seq
		chimp_refs.append(oligo_descriptor)
		# file3.write(oligo_descriptor + '\n')
		# file3.write(seq + '\n')
		pass
		# print x
		# print chimp_ref
		# print modification
	if '-2_modification(s)' in modification:
		seq = x[oligo_descriptor]
		#print oligo_descriptor
		#print seq
		double_modifications.append(oligo_descriptor)
		pass
		# print x
		# print 'single_modification'
		# print modification
	if '-3_modification(s)' in modification:
		seq = x[oligo_descriptor]
		triple_modifications.append(oligo_descriptor)
		pass
print 'hum_refs:', len(hum_refs)
print 'chimp_refs:', len(chimp_refs)
print 'single_modifications:', len(single_modifications)
print 'double_modifications:', len(double_modifications)
print 'triple_modifications:', len(triple_modifications)
# hum_refs: 751
# chimp_refs: 609
# single_modifications: 1978
# double_modifications: 3352

#file2 = open(directory + '513_nonEnhancerFinder_HAR_hum_refs.fa', 'wb')
#file3 = open(directory + '513_nonEnhancerFinder_HAR_chimp_refs.fa' , 'wb')
file4 = open(directory + 'TFBS_mismatch_nonEnhancerFinder_HAR_single_modification.fa', 'wb')
file5 = open(directory + 'TFBS_mismatch_nonEnhancerFinder_HAR_double_modifications.fa', 'wb')
file6 = open(directory + 'TFBS_mismatch_nonEnhancerFinder_HAR_triple_modifications.fa', 'wb')

# for human_ref in hum_refs:
# 	seq = x[human_ref]
# 	print human_ref
# 	print seq
# 	file2.write(human_ref + '\n')
# 	file2.write(seq + '\n')
# file2.close()

# for chimpanzee_ref in chimp_refs:
# 	seq = x[chimpanzee_ref]
# 	print chimpanzee_ref
# 	print seq
# 	file3.write(chimpanzee_ref + '\n')
# 	file3.write(seq + '\n')
# file3.close()

for single_modification in single_modifications:
	seq = x[single_modification]
	print single_modification
	print seq
	file4.write(single_modification + '\n')
	file4.write(seq + '\n')
file4.close()

for double_modification in double_modifications:
	seq = x[double_modification]
	print double_modification
	print seq
	file5.write(double_modification + '\n')
	file5.write(seq + '\n')
file5.close()

for triple_modification in triple_modifications:
	seq = x[triple_modification]
	print triple_modification
	print seq
	file6.write(triple_modification + '\n')
	file6.write(seq + '\n')
file6.close()


########################################################################
ALL FINAL FASTA FILES FOR LIBRARY DESIGN ARE HERE:

# DIRECTORY: /home/hane/FINAL?/FINAL_FINAL/FINAL_FINAL_FINAL

#  5 FILES BEFORE DROPPING oligos less than 50 bps:
# 239_EnhancerFinder_HAR_permutations.fa	--> 4709 oligos
# split_up_CONTROLS.fa --> 253 oligos				    
# TFBS_mismatch_nonEnhancerFinder_HAR_single_modification.fa --> 827 oligos
# MODIFIED_SPLIT_UP_HARs_1067.fa --> 2134 oligos
# TFBS_mismatch_nonEnhancerFinder_HAR_double_modifications.fa -->  1432 oligos
# TFBS_mismatch_nonEnhancerFinder_HAR_triple_modifications.fa --> 1918 oligos

# CREATE 3 REPLICATES of CONTROL and SPLIT_UP_HAR FILES:

# DROP oligos less than 50bp long

# 3X coverage of all types of oligos (except modifications that are not single modifications)
# SEND KATIE TWO FASTA FILES for dropped oligo analysis
# SEND NADAV THE MODIFIED HARS BED FILE