#!/usr/bin/perl -w
=pod
by Alyx Schubert 6/24/10
This script reads in a file with a list of ID's to find. It returns the data
in the file associated with those ID's.
ARGS: The arguments called to the command line are a list of names/IDs file,
	  followed by the file to search for them in.
NOTE: It's easiest when both inputted files are sorted and you know that you
	  will find every ID in the file you are looking for.
NOTE: May need to tailor to whatever kind of file you are parsing through.
	  Current file type set to: *.ngfa
	  
Changes:
7/16/14 find_lines2.pl: Writes over old output instead of adding to existing 
		files with same result name (searchme.txt.subset.txt, eg).  Also added
		the header for the results file from the file searching through.
=cut

my($names, $file, $output, $line1, $line2, $data, $header);
my $ID1 = "NULL";
my $ID2 = "NULL";
my $stat = 0;
my $nomatch = 0; 

while(@ARGV)
{
	$names = shift(@ARGV);
	push(@ARGV, "subset.txt"); #if you want to change the name, change here
	$output = join(".", @ARGV);
	$file = shift(@ARGV);
	shift(@ARGV); #gets rid of the extra name added before
	open (NAMES, "<$names") or die ("Cannot open file $names: $!");
	open (FILE, "<$file") or die ("Cannot open file $file: $!");
	open (OUT, ">$output") or die ("Cannot open file $output: $!");
	
	$header = <FILE>;
	print OUT "$header";
	
	while($line1 = <NAMES>)
	{
		$nomatch = 0;
		$stat = 0;
		
		$line1 =~ /^(\S*)/; #nonwhitespace characters
		$ID1 = $1;		
		print "ID1: $ID1\n";
		print OUT "$ID1\t";
		until($ID1 eq $ID2)
		{
			if(eof(FILE))
			{
				close (FILE) or die ("Cannot close file $file: $!");
				open (FILE, "<$file") or die ("Cannot open file $file: $!");
				$stat++;
				if($stat>1)
				{
					print OUT "Match not found\n";
					$nomatch = 1;
					last;
				}
			}
			$line2 = <FILE>;
			$line2 =~  /(\S*)\t(.*)/; #nonwhitespace characters followed by tab
			$ID2 = $1;
			$data = $2;
			print "ID2: $ID2\n";
		}
		if($nomatch == 1)
		{
			next;
		}
		#$line2 = <FILE>;
		#$line2 =~ /(\S*)\t(.*)/;
		#$data = $1;
		#$data2 = $2;
		#print "\tdata2:$data2\n";
		print OUT "$data\n";

		
		
	} #end while($line1 = <NAMES>)
	
	close (NAMES) or die ("Cannot close file $names: $!");
	close (FILE) or die ("Cannot close file $file: $!");
	close (OUT) or die ("Cannot close file $output: $!");
	
} #end while(@ARGV)