package GTB::Var::VCF;
############################################################

=head1 NAME

GTB::Var::VCF.pm - A Perl module to hold construct VCF-
          formatted text, specifically for bam2mpg genotype
          calls.

=head1 DESCRIPTION

  This module writes out VCF headers and records for MPG
  SNV and DIV calls.

=head1 DATE

 August 1, 2014

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

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

=head1 PUBLIC METHODS

=cut

############################################################
require 5.005_62;
use strict;
use warnings;

use FileHandle;
use Carp;

#our $VERSION  = '$Rev: 7271 $';
our $REVISION = '$Id: VCF.pm 7271 2014-08-31 13:42:48Z nhansen $';

###########################################################

=over 4

=item new()

  This method creates a new GTB::Var::VCF object.  Options set properties 
  for the entire VCF file, whereas properties of individual variant lines
  are passed as options to the "vcf_line" method.

  Input:  -sample_name - name to be written as the header for the sample 
              genotype column.
          -program - program generating the VCF file (default: bam2mpg).
          -fasta_db - GTB::FASTA object for reference sequences (required).
          -shorten - if true, DIV alleles will be shortened, rather than
              widening to unambiguous changes.
          -spec - if true, always write DIV alleles including the reference 
              base one to the left of the inserted or deleted bases (and
              adjust positions accordingly)
          -file_handle - open filehandle for writing (optional).  If none
              is provided, VCF header and record lines are returned as 
              strings.
          -no_ad_tag - if true, suppresses writing of the "AD" format 
              header line.  Use this option if the file will contain no
              allele depth ("AD") information in the sample genotype fields.

  Output: New GTB::Var::VCF object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $sample_name = $params{'-sample_name'} || 'SAMPLE';
    my $program = $params{'-program'} || 'bam2mpg';
    my $fasta_db = $params{'-fasta_db'}
        or die "Must pass GTB::FASTA object for reference as -fasta_db to GTB::Var::VCF constructor!\n";
    my $shorten = $params{'-shorten'};
    my $spec = $params{'-spec'};
    my $file_handle = $params{'-file_handle'};
    my $no_ad_tag = $params{'-no_ad_tag'};

    my $self = { sample_name => $sample_name,
                 program => $program,
                 fasta_db => $fasta_db,
                 shorten => $shorten,
                 spec => $spec,
                 no_ad_tag => $no_ad_tag,
                 file_handle => $file_handle };

    bless $self, $class;

    return $self;

} ## end new

###########################################################

=item vcf_header()

  The method returns a string containing the VCF header, and
  writes it to the file handle for the object, if one was 
  provided.

  Input: GTB::Var::VCF Object
  Output: Scalar string.

=cut

###########################################################
sub vcf_header {
    my $self  = shift;

    my $program = $self->{program};
    my $sample_name = $self->{sample_name};

    my $header_string = '';

    $header_string .= "##fileformat=VCFv4.0\n";
    my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime();
    $mon++;
    $year += 1900;
    $header_string .= sprintf("##fileDate=%d%02d%02d\n", $year, $mon, $mday);
    $header_string .= "##source=$program\n";

    # include info for RM flag for repeat-masked sequence:
    $header_string .= "##INFO=<ID=RM,Number=0,Type=Flag,Description=\"Lower-case reference\">\n";

    # included genotype id's:
    $header_string .= "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    $header_string .= "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    $header_string .= "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth after Quality Filtering\">\n";
    if (!$self->{no_ad_tag}) {
        $header_string .= "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depths for Indexed Alleles after Quality Filtering\">\n";
    }

    # print required eight fields:
    
    $header_string .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";

    if (my $fh = $self->{file_handle}) {
        print $fh $header_string
            or die "Couldn\'t write VCF header string to object's filehandle: $!\n";
    }

    return $header_string;

} ## end reference_base

###########################################################

=item vcf_line()

  The method returns a string containing a VCF record, and
  writes it to the file handle for the object, if one was 
  provided.

  Input: GTB::Var::VCF Object, and the following parameters:
         -chrom	Chromosome for record
         -pos	Position of variant
         -ref_allele	Reference allele for variant
         -genotype	Genotype (alleles separated by ':'
			   for DIV variants).
         -genotype_score	Genotype score (MPG score)
         -variant_score	Variant score (MPV score)
         -coverage	If a reference to an array of align 
                        hashes, each with a "base" value giving 
                        the allele for that alignment, is 
                        passed in this argument, an "AD"
                        tag with the number of reads supporting
                        each allele will be included in the VCF
                        record, as well as a "DP" tag with the
                        total depth.  Otherwise, if a scalar number
                        is passed, no AD tag will be included
                        with the genotype, but the value of
                        "coverage" will be reported with the 
                        "DP" tag.
         -type		Type of variant (SNV or DIV).
  Output: Scalar string.

=cut

###########################################################
sub vcf_line {
    my $self  = shift;
    my %params = @_;

    my $chrom = $params{'-chrom'};
    my $pos = $params{'-pos'};
    my $ref_allele = $params{'-ref_allele'};
    my $genotype = $params{'-genotype'};
    my $genotype_score = $params{'-genotype_score'} || 0;
    my $variant_score = $params{'-variant_score'} || 0;
    my $ra_coverage = (defined($params{'-coverage'})) ? $params{'-coverage'} : 0;
    my $type = $params{'-type'} || 'SNV';

    my $fasta_db = $self->{fasta_db};
    my $shorten = $self->{shorten};
    my $spec = $self->{spec};

    if (!$chrom || !(defined ($pos)) || !(defined($ref_allele)) || !(defined($genotype))) {
        die "Missing required parameter for vcf_line: (-chrom, -pos, -ref_allele, -genotype)\n";
    }

    my $filtered_depth = (ref $ra_coverage eq 'ARRAY') ? @{$ra_coverage} : $ra_coverage;
    my @bases = (ref $ra_coverage eq 'ARRAY') ? map {$_->{'base'}} @{$ra_coverage} : (); 

    # in a repeat (determined by case of reference base)?

    my $uc_ref = uc $ref_allele;
    my $info_flag = ($ref_allele eq $uc_ref) ? '.' : 'RM';
    $ref_allele = $uc_ref;

    my $doc = @bases;
    my $filter = ($doc && $variant_score>=10 && $variant_score/$doc>=0.5) ? 'PASS' : '.';

    my @alleles = ($type eq 'SNV') ? split //, $genotype : split /:/, $genotype;

    my ($lfe, $rfs);
    if ($type eq 'DIV') { # get rid of asterisks, shorten repetitive indels
        $ref_allele = '' if ($ref_allele eq '*');
        foreach my $allele (@alleles, @bases) {
            $allele =~ s:\*::g;
        }
        ($lfe, $rfs) = ($pos =~ /(\d+):(\d+)/) ? ($1, $2) : (undef, undef);
        if (!defined($lfe)) {
            die "Position passed to vcf_line for DIVs must be numbers separated by colon!  Passed value $pos!\n";
        }
        $pos = $lfe + 1;
        # shorten repetitive indels:
        if ( $shorten ) {
            while ( ( $ref_allele ne '' ) && ( !grep { $_ eq '' } @alleles ) ) {
                my $chop = 1; # do we need to chop a base?
                my $last_char = ( $ref_allele =~ /(.)$/ ) ? $1 : '';
                foreach my $allele (@alleles) {
                    my $allele_last_char = ( $allele =~ /(.)$/ ) ? $1 : '';
                    if ( $allele_last_char ne $last_char ) {
                        $chop = 0;
                        last;
                    }
                }
                if ($chop) {
                    chop $ref_allele;
                    foreach my $allele (@alleles, @bases) {
                        chop $allele;
                    }
                }
                else {
                    last;
                }
            }
        }

        # now move coordinate one to the left if necessary:
        if ( ( $ref_allele eq '' ) || ( grep { $_ eq '' } @alleles ) || ( $spec )) {
            $pos--;
            my $ref_base = uc $fasta_db->seq( $chrom, $pos, $pos );
            $ref_allele = $ref_base . $ref_allele;
            foreach my $allele (@alleles, @bases) {
                $allele = $ref_base . $allele;
            }
        }
    }

    my %allele_hash = map { $_ => 1 } @alleles;
    my @nonref_alleles = grep { $_ ne $ref_allele } keys %allele_hash;
    @nonref_alleles = sort @nonref_alleles;
    my $alt_alleles = (@nonref_alleles) ? join ',', @nonref_alleles : '.';
    my @all_alleles = @nonref_alleles;
    unshift @all_alleles, $ref_allele;

    my %index_hash = map { $all_alleles[$_] => $_ } ( 0 .. $#all_alleles );
    my @allele_indices = map { $index_hash{$_} } @alleles;
    @allele_indices = sort {$a <=> $b} @allele_indices;
    my $index_genotype = join '/', @allele_indices;

    my $ad_string;
    if (ref $ra_coverage eq 'ARRAY') {
        my @ad_counts = ();
        for (my $i=0; $i<=$#all_alleles; $i++) {
            my $allele = $all_alleles[$i];
            my $count = grep /^$allele$/i, @bases;
            push @ad_counts, $count;
        }
        $ad_string = join ',', @ad_counts;
    }

    my $format_string = (defined($ad_string)) ? 'GT:GQ:DP:AD' : 'GT:GQ:DP';
    $filtered_depth .= ":$ad_string" if (defined($ad_string));

    my $vcf_line = "$chrom\t$pos\t.\t$ref_allele\t$alt_alleles\t$variant_score\t$filter\t$info_flag\t$format_string\t$index_genotype:$genotype_score:$filtered_depth\n";

    if (my $fh = $self->{file_handle}) {
        print $fh $vcf_line
            or die "Couldn\'t write VCF record string to object's filehandle: $!\n";
    }

    return $vcf_line;

} ## end vcf_line

1;
__END__

=back
