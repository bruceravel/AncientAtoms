#! /usr/bin/perl -w

use strict;

open SGML, $ARGV[0] || die "could not open sgml file $ARGV[0]\n";

my ($base, $pathto, $suffix) = (0,0,0);
require "File/Basename.pm";
($base, $pathto, $suffix) =
  File::Basename::fileparse( $ARGV[0], "sgml" );

open TEXI, ">".$base."texi" || die "could not open $base.texi for output\n";

my %fonts = (
	     "bf" => "\@b",
	     "em" => "\@emph",
	     "it" => "\@i",
	     "sf" => "\@sc",
	     "tt" => "\@code",
	     );
my $fontlabels = join("|", keys(%fonts));

my %chars = (
	     "chi"    => "chi",
	     "delta"  => "delta",
	     "lambda" => "lambda",
	     "lsqb"   => "[",
	     "mu"     => "mu",
	     "num"    => "#",
	     "percnt" => "\%",
	     "pi"     => "pi",
	     "rho"    => "rho",
	     "rsqb"   => "]",
	     );
my $characters = join("|", keys(%chars));

my %lists = (
	     "descrip" => "table",
	     "enum"    => "enumerate",
	     "item"    => "item\n",
	     "itemize" => "itemize",
	     "quote"   => "quotation",
	     );
my($alllists) = join("|", keys(%lists));


while (<SGML>) {

    chomp;
    (/<tag\/(.+)\//) && do {
	my($replacement) = "\@item $1";
	s/$&/$replacement/;
    };
    my(@line) = split;
    my($word) = 0;
    foreach $word (@line) {
	MATCH : {
	                               # sgml-ized special characters
	    ($word =~ /&($characters);/) && do {
		my($replacement) = $chars{$1};
		$word =~ s/$&/$replacement/;
		last MATCH;
	    };
	                               # replace <p> with nothing
	    ($word =~ /(<p>)/) && do {
		$word =~ s/$1//;
		last MATCH;
	    };
	                               # <tscreen><verb> 
	    ($word =~ /(<tscreen><verb>)/) && do {
		$word =~ s/$1/\@example/;
		last MATCH;
	    };
	    ($word =~ /(<\/verb><\/tscreen>)/) && do {
		$word =~ s/$1/\@end example/;
		last MATCH;
	    };

	    
	                               # list structures
	    ($word =~ /<(\/?)($alllists)>/) && do {
		my($replacement) = $1 ? "\@end " : "\@";
		$replacement .= $lists{$2};
		$word =~ s/$&/$replacement/;
		last MATCH;
	    };

	                                # fonts
	    ($word =~ /<($fontlabels)\/(.+)\//) && do {
		my($replacement) = $fonts{$1} . "\{" . $2 . "\}";
		$word =~ s/$&/$replacement/;
		last MATCH;
	    };

	                                # sectioning
	    ($word =~ /(<sect>)/) && do {
		$word =~ s/$1/\@chapter/;
		last MATCH;
	    };
	    ($word =~ /(<sect1>)/) && do {
		$word =~ s/$1/\@section/;
		last MATCH;
	    };
	    ($word =~ /(<sect2>)/) && do {
		$word =~ s/$1/\@subsection/;
		last MATCH;
	    };

	}
	print TEXI "$word ";
    }
    print TEXI "\n";
}


