#!/usr/bin/perl
#
# mkplot.pl
# takes the strange output from funpoints and creates a file suitable for
# plotting with gnuplot

# catch ctrl-c
$SIG{INT} = sub { die "\nSignal caught, exiting.\n" };

# open pipe to gnuplot
open(GNUPLOT,"|gnuplot") || die "Could not open pipe to gnuplot.\n";
# set initial plot parameters
print GNUPLOT "set terminal X11\n";
print GNUPLOT "unset key\n";
print GNUPLOT "set size ratio 1\n";
print GNUPLOT "set xrange [0:1]\n";
print GNUPLOT "set yrange [0:1]\n";

# loop forever, until ctrl-c is caught
for(;;){
	# ptfifo is the named pipe you cat the output of funpoints to
	open(FIFO,"cat ptfifo |") || die "Could not open named input.\n";
	@a = <FIFO>;
	# holds all the lines
	my @lines;

	# go through the contents...
	foreach(@a){
		# lines beginning with hashes are comments and are ignored
		if(m/^.*#/){ next; }

		# convert all non-numeric fluff to white space
		$_ =~ s/(:|\)|\(|-|>)//g;

		# no need for that newline
		chomp;
		# split based on whitespace and append to other array
		@tmp= split;
		push @lines, [ @tmp ];
	}

	# gnuplot accepts '-' as a plot option, assuming the data will follow one
	# line at a time.  The EOF for a set of data is the letter 'e' on a line by
	# itself.  So for N line segments, make N seperate plots and choose the
	# linestyle (option 'ls') based on the grouping.  This assures the colors
	# will be consistent.
	print GNUPLOT "plot '-' w lp ls $lines[0]->[7]";
	for(my $i=1;$i<scalar(@lines);++$i){
		print GNUPLOT ", '-' w lp ls $lines[$i]->[7]";
	}
	print GNUPLOT "\n";

	# now print a the list of lines
	foreach(@lines){
		print GNUPLOT "$_->[2]\t$_->[3]\n";
		print GNUPLOT "$_->[4]\t$_->[5]\n";
		print GNUPLOT "e\n";
	}
	print GNUPLOT "\n";
	# close up and do it again
	close(FIFO);
}
