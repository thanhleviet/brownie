#!/usr/bin/perl -w
use diagnostics;
open (IN,"cdfvectorholder.cpp");
open (OUT,">cdfvectorholderMOD.cpp");
$repeatcount=0;
$repeatphrase="";
while (<IN>) {
	$inline=$_;
	chomp $inline;
	if ($inline=~m/(contentsofrow.push_back\(1e\-05\)\; \/\/obs value was actually 0)/) {
		$repeatcount++;
		$repeatphrase=$1;
	}
	elsif ($inline=~m/(contentsofrow.push_back\(1\)\;)/) {
		$repeatcount++;
		$repeatphrase=$1;
	}
	else {
		if ($repeatcount>1) {
			print OUT "    for (int i=0;i<$repeatcount;i++) {\n        $repeatphrase\n    }\n";
			$repeatcount=0;
		}
		elsif ($repeatcount==1) {
			print OUT "    $repeatphrase\n";
			$repeatcount=0;
		}
		print OUT "$inline\n";
	}
}
close IN;
close OUT;