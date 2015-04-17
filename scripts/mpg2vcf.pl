#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
# $Id: mpg2vcf.pl 7267 2014-08-29 20:59:27Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(:ALL);
use GTB::FASTA;
use GTB::Var::VCF;
use GTB::Unix qw(which);
use GTB::File qw(Open :use_bgzip);
use vars qw($VERSION);

our %Opt;

my $Id = q$Id: mpg2vcf.pl 7267 2014-08-29 20:59:27Z nhansen $;
$VERSION = sprintf "%.4f", substr( q$Revision: 0$, 10 ) / 10000;

=head1 NAME

mpg2vcf.pl - script to convert bam2mpg output to VCFv4.0

=head1 SYNOPSIS

  mpg2vcf.pl --mpg <mpg.out> --ref <fasta> --sample <sample_name> <--shorten> [-notabix]

or

  mpg2vcf.pl --mpgfof <file listing mpg files> --ref <fasta> --sample <sample name> <--shorten> [-notabix]

=head1 DESCRIPTION

This script reads in MPG output from bam2mpg and writes a VCFv4.0 file containing the same calls.  If output files end in '.gz' it will compress with bgzip and create tabix indexes. To prevent index creation pass -notabix option

=cut

#------------
# Begin MAIN
#------------

process_commandline();
my @mpg_files = split("," , $Opt{'mpg'} );
my $ref_fasta   = $Opt{'ref'};
my $sample_name = $Opt{'sample'} || 'Unknown';
my $snv_outfile = $Opt{'snv_outfile'};
my $div_outfile = $Opt{'div_outfile'};

my $samtools_exe = 'samtools';
$GTB::FASTA::SAMTOOLS = $samtools_exe;

my $fasta_db = GTB::FASTA->new($ref_fasta);

my $mode = $Opt{'append'} ? 'a' : 'w';
my $snv_fh = Open($snv_outfile, $mode);
my $snv_vcf_obj = GTB::Var::VCF->new(-sample_name => $sample_name,
                                     -program => 'mpg2vcf.pl',
                                     -fasta_db => $fasta_db,
                                     -shorten => $Opt{'shorten'},
                                     -spec => $Opt{'spec'},
                                     -no_ad_tag => 1,
                                     -file_handle => $snv_fh);


my ($div_fh, $div_vcf_obj);
if ( $div_outfile eq $snv_outfile ) {
	$div_fh = $snv_fh;
        $div_vcf_obj = $snv_vcf_obj;
}
else {
	$div_fh = Open($div_outfile, $mode);
        $div_vcf_obj = GTB::Var::VCF->new(-sample_name => $sample_name,
                                          -program => 'mpg2vcf.pl',
                                          -fasta_db => $fasta_db,
                                          -shorten => $Opt{'shorten'},
                                          -no_ad_tag => 1,
                                          -file_handle => $div_fh);
}

$snv_vcf_obj->vcf_header() unless ( $Opt{'noheader'} );
$div_vcf_obj->vcf_header() if ( ( $snv_outfile ne $div_outfile )
	                             && ( !$Opt{'noheader'} ) );

for my $file (@mpg_files){
  my $mpg_fh = Open($file);

  while (<$mpg_fh>) {
	if ( /^MPG(?:_SNV)?\s(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\d+)\s([01])
         (?:\s(\d+))?(?:\s(\d+))?$/x) {
		my ( $chr, $pos, $ref_allele, $genotype, $mpg, $ref_nonref,
			$coverage, $mpv ) = ( $1, $2, $3, $4, $5, $6, $7, $8 );
                my $score = defined($mpv) ? $mpv : (($ref_nonref) ? $mpg : 0);
		if ($Opt{'min_score'} && $score < $Opt{'min_score'}) {
			next;
		}
                if ($Opt{'variants'} && $ref_nonref == 0) {
                        next;
                }
                my $vcf_line = $snv_vcf_obj->vcf_line(  '-chrom' => $chr,
                                                        '-pos' => $pos,
                                                        '-ref_allele' => $ref_allele,
                                                        '-genotype' => $genotype,
                                                        '-genotype_score' => $mpg,
                                                        '-variant_score' => $score,
                                                        '-coverage' => $coverage,
                                                        '-type' => 'SNV' );
	}
	elsif ( /^MPG(?:_DIV)?\s(\S+)\s(\d+):(\d+)\s(\S+)\s(\S+)\s(\d+)\s([01])
             (?:\s(\d+))?(?:\s(\d+))?$/x) {
		my ( $chr, $lfe, $rfs, $ref_allele, $genotype, $mpg, $ref_nonref,
			$coverage, $mpv )
		  = ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
                my $score = defined($mpv) ? $mpv : (($ref_nonref) ? $mpg : 0);
		if ($Opt{'min_score'} && $score < $Opt{'min_score'}) {
			next;
		}
		if ( ( !$Opt{'keep_ref_div'} ) && ( $ref_nonref == 0 ) ) {
			next;
		}
		next if ( $lfe < 1 );    # something weird going on!
                my $vcf_line = $div_vcf_obj->vcf_line(  '-chrom' => $chr,
                                                        '-pos' => "$lfe:$rfs",
                                                        '-ref_allele' => $ref_allele,
                                                        '-genotype' => $genotype,
                                                        '-genotype_score' => $mpg,
                                                        '-variant_score' => $score,
                                                        '-coverage' => $coverage,
                                                        '-type' => 'DIV' );
	}
    }
}
close($snv_fh) or die "Error closing $snv_outfile, $!\n";
close($div_fh) or die "Error closing $div_outfile, $!\n";

if ( $snv_outfile =~ /gz$/ && $Opt{tabix} ) {
    unlink "$snv_outfile.tbi";  # tabix won't overwrite, so delete first
	my $cmd = " $Opt{tabix} -p vcf $snv_outfile ";
	my $w = `$cmd`;
}
if ( $div_outfile =~ /gz$/ && $Opt{tabix} ) {
    unlink "$div_outfile.tbi";
	my $cmd = " $Opt{tabix} -p vcf $div_outfile ";
	my $e = `$cmd`;
}

#------------
# End MAIN
#------------

sub process_commandline {
	# Set defaults here
	%Opt = ( ref => '/scratch/fasta/hg18/hg18.mfa', );
	GetOptions(
		\%Opt, qw(
		  manual help+ version
		  verbose ref=s mpg=s
		  mpgfof=s keep_ref_div variants
		  sample=s snv_outfile=s
		  div_outfile=s append_char gzip
		  noheader shorten tabix! min_score=i
                  spec
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "mpg2vcf.pl, ", q$Revision: 7267 $, "\n";
	}

	if ( $Opt{mpgfof} ) {
		my $fh = Open($Opt{mpgfof});
		my @mpg_files = ();
		while (<$fh>) {
			chomp;
			push @mpg_files, $_;
		}
		$Opt{mpg} = join ',', @mpg_files;
	}
	if ( !exists( $Opt{mpg} ) ) {
		die "Specify input mpg filename with option --mpg.\n";
	}
	if ( !$Opt{sample} ) {
		my ( $volume, $path, $file ) = splitpath( $Opt{mpg} );
		if ( $file =~ /^(.*?)\./ ) {
			$Opt{sample} = $1;
		}
		else {
			$Opt{sample} =~ s/\.gz//;    # in case it is not mpg.out
		}
	}

	if ( !$Opt{snv_outfile} ) {
		$Opt{snv_outfile} = "$Opt{sample}.mpg.snv.vcf";
		if ($Opt{gzip}) {
			$Opt{snv_outfile} .= ".gz";
		}
	}

	if ( !$Opt{div_outfile} ) {
		$Opt{div_outfile} = "$Opt{sample}.mpg.div.vcf";
		if ($Opt{gzip}) {
			$Opt{div_outfile} .= ".gz";
		}
	}
	if ( !defined $Opt{tabix} ) {
		$Opt{tabix} = which("tabix");
		chomp $Opt{tabix};
	}
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--mpg> file.mpg.out

Specify the MPG output file to be converted to VCFv4.0.

=item B<--mpgfof> file_of_files

Specify a file that lists all of the MPG output files to use as input.

=item B<--ref>

Specify the reference fasta file that was used in the bam2mpg run that resulted
in the MPG output file.

=item B<--sample>

Specify the sample name from which these genotypes come.  This will be used
inside the VCF file as the column header, and as the base name for the
default output file names.

=item B<--snv_outfile>

Specifies the name of the output VCF file for SNV variants.  If the filename ends
in "gz", will zip with bgzip. (default: sample.mpg.snv.vcf)

=item B<--div_outfile>

Specifies the name of the output VCF file for DIV variants.  If the filename ends
in "gz", will zip with bgzip.  If it is the same as the file specified in
--snv_outfile, will write all variants to a single file.
(default: sample.mpg.div.vcf)

=item B<--variants>

When the --variants option is specified, VCF lines will only be output for variants.
(default: print lines for all sites, even those with homozygous reference genotypes.)

=item B<--shorten>

When the --shorten option is specified, repetitive DIV variants will be shortened
to their shortest variant. (Note: currently, vcf2mpg.pl will not re-create the
original MPG file when this option is used).

=item B<--notabix>

By default when a output file ends in 'gz', tabix will create the corresponding index.
To override this behaviour specify --notabix.

=item B<--min_score> N

Minimum MPG score to create VCF record.  Variants with scores below this
threshold will not be included.  Default is no minimum score.

=item B<--gzip>

Generate compressed output files (suffixed with .gz extension).  Has no
effect when C<--snv_outfile> or C<--div_outfile> are explicitly specified.

=item B<--spec>

Print indel lines exactly according to the V4.1 VCF specifications.  That is, position
the insertion or deletion (or mixed type variant) one base to the left of the "variant",
and thus including the reference base in that position in both the reference allele and
all alternate alleles.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use.

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose.

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation.

=cut
