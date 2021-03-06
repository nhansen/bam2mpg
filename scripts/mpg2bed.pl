#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
# $Id: mpg2bed.pl 7202 2014-07-25 12:35:04Z nhansen $

use warnings;
use strict;
use GTB::File qw(Open);
use Getopt::Long;
use Pod::Usage;
our (%Opt, $OutFh);

my @in_files;
my $last_chrom = "";
my $last_pos = "";
my $start_pos = "";

process_commandline();
for my $file (@ARGV) {
    if ($file =~ /\.fof$/) {
        my $ffh = Open($file);
        while (<$ffh>) {
            chomp;
            push @in_files, $_;
        }
    }
    else {
        push @in_files, $file;
    }
}

$OutFh = Open($Opt{output}, 'w');
print $OutFh "track name=\"$Opt{sample}\" description=\"$Opt{sample} - " .
    ($Opt{mpv} ? "Bases interrogated (MPV score >= $Opt{min_score} or MPG score >= $Opt{min_score}"
               : "Bases confidently called (MPG score >= $Opt{min_score}");
if ($Opt{ratio}) {
    print $OutFh " and score/depth ratio >= $Opt{ratio}";
}
if ($Opt{depth}) {
    print $OutFh " and depth of coverage >= $Opt{depth} reads";
}
print $OutFh qq{)"\n};

for my $file (@in_files) {
    my $fh = Open($file);
    while (<$fh>) {
        if (/^MPG(?:_SNV)?\s(\w+)\s(\d+)\s\w\s\w+\s(\d+)\s([01])\s(\d+)(?:\s(\d+))?/) {
            my ($chrom, $pos, $score, $nonref, $covg) = ($1, $2, $3, $4, $5);
            next if $Opt{nonref} && !$nonref;
            if ($Opt{mpv} && $6) {  # when MPV score > 0, it is >= MPG score
                $score = $6;
            }
            my $ratio = $covg > 0 ? $score / $covg : 0;
            next if ($score < $Opt{min_score} || $ratio < $Opt{ratio} );
            next if ($Opt{depth} && $covg < $Opt{depth} );
            if ($chrom ne $last_chrom || ($pos != $last_pos + 1) ) {
                if ($last_chrom) {
                    print $OutFh join("\t", $last_chrom, $start_pos - 1,
                            $last_pos), "\n";
                }
                $start_pos = $pos;
                $last_chrom = $chrom;
            }
            $last_pos = $pos;
        }
    }
}
if ($last_chrom) {
    print $OutFh join("\t", $last_chrom, $start_pos - 1, $last_pos), "\n";
}

sub process_commandline {
    #Defaults
    %Opt = (
         min_score => 10,
         output  => '-',
         ratio   => 0,
         depth   => 0,
    );
    GetOptions(\%Opt, qw(
        depth|mincovg=i
        help
        manual
        min_score|threshold|score=i
        mpg_file=s
        mpv
        nonref
        output=s
        ratio=f
        sample|name=s
    )) || pod2usage();
    if ($Opt{manual}) { pod2usage(verbose => 2); }
    if ($Opt{help})   { pod2usage(verbose => 1); }
    if ($Opt{mpg_file}) {
        unshift @ARGV, $Opt{mpg_file};
    }
    if (!@ARGV) {
        pod2usage("At least one MPG file is required");
    }
    if (!$Opt{sample}) {
        $Opt{sample} = $ARGV[0];
    }
}

=head1 NAME

mpg2bed.pl - Write BED file to represent coverage of bam2mpg calling

=head1 SYNOPSIS

Report regions with a confident (MPG >= 10, MPG/depth ratio >= 0.5) genotype:

    mpg2bed.pl -ratio 0.5 -sample NAME -o NAME.mpg.bed.gz *.mpg.out

Report regions that have sufficient coverage for MPG to determine whether
there is a variant present at the position (MPV >=10 or MPG >=10):

    mpg2bed.pl -mpv -sample NAME -o NAME.mpv.bed.gz *.mpg.out

See C<mpg2bed.pl -h> for all options.

=head1 DESCRIPTION

This script converts MPG output to a BED file representation.  It includes in
the bedfile positions at which the MPG (or MPV) score is greater than or
equal to a threshold.  Adjacent positions are grouped into the same bed line.

=head1 INPUT

One or more MPG output files may be provided as arguments.  Or one file may
be specified with the C<--mpg_file> option.

=over 4

=item B<--mpg_file> FILE

The file provided should be an MPG output file, generated by bam2mpg.  This
file may also be a "file of files" (must have a .fof extension) containing
absolute paths to a set of MPG output files, i.e. MPG output broken up by
chromosome.

=back

=head1 OPTIONS

=over 4

=item B<--depth> NUMBER

Minimum depth of coverage, as reported by bam2mpg (typically, reads with Q20
bases covering the position, capped at 200).

=item B<--min_score> SCORE

Minimum score (MPG or MPV) value for a position to be considered "covered".

=item B<--mpv>

Use MPV score, rather than MPG genotype score to determine which bases are
covered; C<--min_score> (and C<--ratio>, if provided) will apply to the MPV
score instead.

=item B<--nonref>

Output only positions where the most probable genotype is not homozygous for
the reference allele.

=item B<--output> /path/to/output.bed[.gz]

Specify path for output file (otherwise, will write to STDOUT).

=item B<--ratio> NUMBER

The minimum score / coverage ratio required to consider a position
"covered", and thus be included in the output BED file.

=item B<--sample> NAME

A name for the sample, which will be written to the "track" line.

=back

=head1 AUTHORS

 Jamie K. Teer <teerj@mail.nih.gov>
 Peter Chines <pchines@mail.nih.gov>

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
authors as the source of the software or data should be made.

=cut
