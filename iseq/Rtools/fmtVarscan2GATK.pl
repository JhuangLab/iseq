#!perl

##Use third party modules
use diagnostics;
use strict;
use warnings;
use Getopt::Long;

## check threads available or not
$| = 1;

my $input;
my $output;

GetOptions(
'i|input|=s' => \$input,
'o|output|=s' => \$output,
);


my $output_file = "${output}.tmp";
open RES,">$output_file" || die "Can\'t output to file $output_file:$!\n";

my $filename = $input;

open (T,"$filename") || die "Can\'t open $filename:$!\n";

while (<T>) {
	chomp;
	if (/^#/){
		print RES "$_\n";
	}else{
		my @n = split "\t",$_;
		my $n5 = $n[4];

		my $out1;
		my $out2;
		if ($n5 =~ /(.*?)\/(.*)/){
			my $d1 = $1;
			my $d2 = $2;

			if ($d1 =~ /[+-]/){
				#nothing
			}else{
				splice @n,4,1,$d1;
				$out1 = join "\t",@n;
				print RES "$out1\n";
			}
			if ($d2 =~ /[+-]/){
				#nothing
			}else{
				splice @n,4,1,$d2;
				$out2 = join "\t",@n;
				print RES "$out2\n";
			}
		}else{
			print RES "$_\n";
		}

	}

}

close (T);
close RES;

$output_file = $output;
open RES,">$output_file" || die "Can\'t output to file $output_file:$!\n";

$filename = "${output}.tmp";

open (T,"$filename") || die "Can\'t open $filename:$!\n";

while (<T>) {
	chomp;
	if (/^#/){
		print RES "$_\n";
	}else{
		my @n = split "\t",$_;
		my $n4 = $n[3];

		my $out1;
		my $out2;
		if ($n4 =~ /(.*?)\/(.*)/){
			my $d1 = $1;
			my $d2 = $2;

			if ($d1 =~ /[+-]/){
				next;
			}else{
				splice @n,3,1,$d1;
				$out1 = join "\t",@n;
				
			}
			if ($d2 =~ /[+-]/){
				#nothing
			}else{
				splice @n,3,1,$d2;
				$out2 = join "\t",@n;
				
			}
			defined $out1 ? print RES "$out1\n" : print RES "$out2\n";
		}else{
			print RES "$_\n";
		}

	}

}

close (T);
close RES;

system "rm ${output}.tmp";
