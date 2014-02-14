#!/usr/bin/env perl

### Written by Kevin Weitemier, copyright(c) 2011, version 0.0, 17 Feb 2011
### This program is free software: you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free Software
### Foundation, either version 3 of the License, or (at your option) any later
### version. A copy of this license is available at <http://www.gnu.org/licenses/>.
###	version 0.01 23 Feb 2011 -	made it slightly faster, maybe?
###	version 1.00b 03 Mar 2011 -	Added ability to handle lines where the consensus differs from the reference.
###	version 1.01b 09 Mar 2011 -	Fixed a bug when the sequence line contains only Ns.
###	version 2.00b 11 Apr 2011 -	Reduced memory use, added output fields, fixed bug related to ^$.
###	version 3.00b 14 Apr 2011 -	Added section displaying the differing reads and their qual scores.
###	version 3.01b 15 Apr 2011 -	Explicitly undefined variables, which maybe helps memory usage.
###	version 3.02b 18 Apr 2011 - 	Prints a shortened but otherwise unmodified reads line.
###	version 3.03b 29 Oct 2013 -	Modified to work with mpileup output
### Purpose: Calculates the proportion of reads at a position
### 	that differ from the majority.
### Syntax: polymorphic_read_counter_bwaPileup.pl [filename]
### [filename] must be a file in samtools mpileup format.
### Output is tab-delimited file with:
### Name, pileup position, reference base, proportion of non-consensus reads, number of reads, consensus base, differing bases (may be blank), their qual scores (may be blank), all bases, all qual scores.

use diagnostics;
use strict;

while (<>) {					### This goes through the input file, line by line.
	chomp;
	my @fields = split(/\t/,$_);			### Splits the line into fields based on tabs, puts them in @fields
	unless($fields[2] eq "*") {			### Unless the third field contains a * (is an indel line), do this operation
		my $line = $fields[4];			### Assigns the reads field to $line
		my $differ_count;			### Initializes this variable
		$line =~ s/\^.//g;			### Removes all ^ signs and the following character
		$line =~ s/\$//g;			### Removes all $ signs

		while ($line =~ m/-(\d+)/) {		### If the line contains a - sign then a number (including 10 or more), performs operation and stores the number as $1
			$line =~ s/-$1.{$1}//;		### Removes the -digit and the following # of characters equal to the digit.
		}

		while ($line =~ m/\+(\d+)/) {		### Performs the same operation but with + then a number
			$line =~ s/\+$1.{$1}//;		
			#$differ_count += 1;		### Adds the number of reads with insertions to $differ_count.# Realized this can't really be used because insertions can't be aligned to reference. Still, could be interesting.
		}

		my $short_line = $line;			### Keeps the shortened line to be printed later, so not affected by later conversions for calculations (i.e. a -> A).
		my $length = length($line);				### Calculates the length of this variable (the read depth).
		$differ_count += ($line =~ tr/ACTGactg*/ACTGACTG*/);	### Calculates the number of positions with those characters
		my $prop = ($differ_count / $length);			### Calculates the proportion of characters that are different
###################### This section displays bases differing from the reference and their corresponding quality scores #######################################
		my $qual_length = length($fields[5]);			### Finding the length of the qual line
		my $differ_reads = "";					### Initiating variables
		my $differ_quals = "";
		my $different_bases_info = "";

		if ($length == $qual_length) {				### Testing if the read and qual lines are the same length.
			while ($line =~ m/[^,\.]/g){			### Matches characters that aren't , or .
                	        $differ_reads .= $&;			### Appends these characters to this variable.
        	                $differ_quals .= substr($fields[5], (pos($line)-1), 1);		### Appends the corresponding character from $fields[5] to $differ_quals
	                }
			$different_bases_info = $differ_reads . "\t" . $differ_quals;		### Puts these together, to be printed later
		} else { $different_bases_info = "These fields are NOT equal!!! (Meaning something is probably wrong with the line input and modification.)";}
###################### It's possible for a non-reference base to be the most numerous, but less than 50% due to N's, so this next section handles that. #################
		my $not_same_count = ($line =~ tr/ACTG*Nn/ACTG*Nn/);
		my $Nprop = ($not_same_count / $length);		### This and the previous line calculate the proportion of reads not matching the reference (which includes N's).

		if ($Nprop >= 0.5) {					### If $Nprop > 0.5, the consensus could differ from the reference, requiring different calculations (but would take longer to perform on all lines).
			$line =~ s/\./,/g;				### Converts dots . to commas , so they can be counted together.
			my @char_array = split(//,$line);		### Splits the sequence into individual characters	
			my %char_hash = ( 'n' => 0 );				### Initiates this hash
			my $consensus = 'n';					### Initiates this variable
			foreach my $char (@char_array) {			### This loop sets $consensus equal to the most numerous character in the sequence...
				unless ($char =~ /[nN]/) {			### ... unless it's an N or n
					$char_hash{$char} += 1;		### This adds up each character occurance (by adding 1 to the hash value for that key).
					$consensus = $char if ($char_hash{$char} > $char_hash{$consensus}) 	### If this character has been seen more times than the current consensus character make $char the new $consensus.
				}
			}
			my $n_count = ($line =~ tr/Nn/Nn/);					### Counts all N's.
			my $new_prop = (($length - $char_hash{$consensus} - $n_count)/$length);	### Calculates the proportion of differing reads.
			print $fields[0],"\t",$fields[1],"\t",$fields[2],"\t",$new_prop,"\t",$length,"\t",$consensus,"\t",$different_bases_info,"\t",$short_line,"\t",$fields[5],"\n";
			undef @char_array; undef %char_hash; undef $consensus; undef $n_count; undef $new_prop;
		} else {
			print $fields[0],"\t",$fields[1],"\t",$fields[2],"\t",$prop,"\t",$length,"\t",$fields[2],"\t",$different_bases_info,"\t",$short_line,"\t",$fields[5],"\n";
		}
		undef $line; undef $differ_count; undef $short_line; undef $length; undef $prop; undef $qual_length; undef $differ_reads; undef $differ_quals; undef $different_bases_info; undef $not_same_count; undef $Nprop;
	}
}
### EOF ###
