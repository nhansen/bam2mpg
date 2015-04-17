# $Id:$
package GTB::Unix;

use strict;
use File::Temp qw(tempfile);
use Carp qw(carp croak confess);
use FileHandle;
use File::Basename;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( which tailfile grepstring cp );

sub which {
    my ($prog) = @_;

    my @p = split ':', $ENV{PATH};
    for my $p (@p) {
        if (-f "$p/$prog" && -x "$p/$prog") {
            return "$p/$prog";
        }
    }
    return;
}

# perl implementation of tail command, based on code at:
# http://www.thbz.org/ppt/tail

sub tailfile {
    my $file = shift @_;
    my $no_lines = shift @_;

    if ($no_lines !~ /^\d+$/) {
        die "Second argument to tailfile must be a number of lines!\n";
    }

    my $fh = new FileHandle("$file", "r")
        or die "Couldn\'t open $file: $!\n";

    # Start from the end, read the $pos last characters, count
    # how many lines it makes, read previous $pos characters
    # if necessary, and so on, until -$no_lines lines have been read

    my @tail; # will store final $no_lines lines

    my $nl=0;

    seek $fh, 0, 2; # seek to end of file

    my $step = tell $fh; # current pos
    my $pos = -8;
    my $tmp;
    my @slice;
    my $last = 0;
	  
    REWIND: while (1) {
        $pos *= 2;
 
        unless(seek $fh, $pos, 2) {
            seek $fh, 0, 0;
            $last = 1; # hit beginning of file
        } else {
            next REWIND unless <$fh>; # Ignore the current line
        }
 
        $tmp = tell $fh;
        @slice = ();
	     
        READLINE: while (<$fh>) { # read up to the previous spot we were at, pushing lines onto slice
            $nl++;
	    push @slice, $_;
	    last READLINE if (tell $fh >= $step);
        }
	     
        push @tail, reverse @slice;
	     
        last REWIND if ($nl >= $no_lines or $last == 1);
        $step = $tmp;
    }

    close $fh;

    my @return_list = ();
    while (@return_list < $no_lines && @tail) {
        unshift @return_list, shift @tail;
    }

    return @return_list;
}

sub grepstring {
    my $file = shift @_;
    my $string = shift @_;

    my $fh = new FileHandle("$file", "r")
          or die "Couldn\'t open $file: $!\n";
    my @return_lines = grep /$string/, <$fh>;

    return @return_lines;
}

sub cp {
    my $source = shift;
    my $dest   = shift;
    die "cp() cannot yet copy multiple files" if (@_);
    die "cp() cannot copy directories" if (-d $source);
    open(SRC, "$source")
        or die "Couldn't open $source: $!";
    if (-d $dest) {
        my $filename = basename($source);
        $dest .= "/$filename";
    }
    open(DEST, ">$dest") or die "Couldn't O_CREAT $dest: $!";
    my $buf;
    while(read(SRC, $buf, 4096)) {
        print DEST $buf;
    }
    close(SRC);
    close(DEST);
}

1;
__END__

=head1 NAME

GTB::Unix - Functions to replace Unix shell commands

=head1 SYNOPSIS

    use GTB::Unix qw(which);
    my $path_to_program = which($program);

=head1 DESCRIPTION

This is a library of exportable functions.

=head1 EXPORTABLE FUNCTIONS

=head2 Unix commands:

=head3 which($program)

The which command traverses the directories in the "PATH" environment 
variable (delimited by the ":" character), and returns the first 
instance of an existing executable file of the form $dir/$program.

  Usage : my $fullpath = which('command');
  Args  : command to be searched for in user's path
  Return: scalar string with full path of program

=head3 tailfile($file, $no_lines)

Returns the last $no_lines lines of $file in an array (in the order they
exist in the file, with the last line of the file in the highest indexed
element of the array).  As with the Unix tail command, if there are fewer
lines in the file than the desired number of lines, will only return as
many lines as are in the file.

  Usage : my @last_lines = tailfile('myfile.txt');
  Args  : file name, number of desired lines (which must be an integer)
  Return: array of file lines
  Except: dies if file does not exist or second argument is not an integer.

=head3 grepstring($file, $mystring)

Returns an array of lines from $file containing the string $mystring.  Note
that this function does NOT consider regular expression.

  Usage : my @last_lines = grepstring($file, $mystring);
  Args  : file name, string to search for
  Return: array of file lines (eventually would be nice to return a single
          concatenated string if called in scalar context).
  Except: dies if file does not exist

=head3 cp($from, $to)

Copies a file from one place to another. Dies on errors. This function does
not handle copying multiple files to a directory the way the system command does.

  Usage : cp($from, $to);
  Args  : source file, destination file
  Return: None
  Except: Dies if file cannot be copied for some reason.

=head1 AUTHOR

 Peter Chines <pchines@mail.nih.gov>
 Nancy Hansen <nhansen@mail.nih.gov>

=head1 LEGAL

This software is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S. Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
