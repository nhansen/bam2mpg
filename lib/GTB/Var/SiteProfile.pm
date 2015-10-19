package GTB::Var::SiteProfile;
############################################################

=head1 NAME

GTB::Var::SiteProfile.pm - A Perl module to hold information about 
          bases of short sequencing reads aligned to a
          particular position of a reference sequence.

=head1 DESCRIPTION

  This module allows the user to store, retrieve, and 
  manipulate the data from short sequencing reads at a 
  particular position in a reference sequence.

=head1 DATE

 February 8, 2007

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

our $VERSION  = '0.01';
our $REVISION = '$Id: SiteProfile.pm 7230 2014-08-02 00:49:21Z nhansen $';

###########################################################

=over 4

=item new()

  This method creates a new GTB::Var::SiteProfile object.

  Input:  -reference_base (or -reference_allele for DIPs) - the 
              base or bases within the reference sequence at 
              this position.
          -coverage - a reference to an array of bases at
              this position, where each element is a reference
              to a hash with the following fields:
                  base => $base (base identity, required)
                  error_prob => $prob (probability this base 
                         has been read incorrectly. Optional)
                  solexa_scores => $ra_scores (reference to
                         an array of four Solexa-style scores,
                         for bases A,C,G, and T respectively.
                         Optional)
          -reference - name of the reference sequence (Optional).
          -position - position within the reference sequence 
              (1-based. Optional)
          -default_error - value of the "error_prob" to use if
                  neither error_prob nor solexa_scores is 
                  populated in a coverage element (default 0.01).

          FOR HANDLING DIPS:

          -allele_seqs - reference to an array of all possible
              allele sequences, one of which can be an empty 
              sequence ('').  When this parameter is passed to
              the constructor, it overrides the default, single-base
              genotypes (AA, AC, etc.) with all possible combinations
              of the alleles passed.
          
          FOR HANDLING SINGLE COPY REGIONS:

          -single_copy - if true, will calculate genotypes that
              consist of a single allele, e.g., A, T, G, or C
              for SNVs.
          -verbose - if true, will turn on various verbose reporting.

  Output: New GTB::Var::SiteProfile object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $reference_allele = defined ($params{-reference_base}) ?
                           $params{-reference_base} : $params{-reference_allele};
    my $ra_coverage = $params{-coverage} || [];
    my $reference = $params{-reference};
    my $position = $params{-position};
    my $default_error = $params{-default_error} || 0.01;
    my $ra_allele_seqs = $params{-allele_seqs};
    my $single_copy = $params{-single_copy};
    my $verbose = $params{-verbose};

    my @genotypes = ($single_copy) ? qw( A C G T ) :
                               qw( AA AC AG AT CC CG CT GG GT TT );

    if ($ra_allele_seqs) # user has specified alleles--create new genotypes
    {
        # order alleles:
        my @ordered_alleles = sort @{$ra_allele_seqs};
        $ra_allele_seqs = \@ordered_alleles;
        @genotypes = ();

        if ($single_copy) {
            @genotypes = @ordered_alleles; 
        }
        else {
            for (my $i = 0; $i<= $#{$ra_allele_seqs}; $i++) # first allele
            {
                my $allele1 = $ra_allele_seqs->[$i];
                if (($allele1 ne '') && ($allele1 ne '*') && ($allele1 eq lc $allele1))
                {
                    carp "You\'ve passed a lower case allele $allele1 to SiteProfile.  Is this what you meant to do?\n";
                }
                for (my $j = $i; $j<= $#{$ra_allele_seqs}; $j++) # second allele
                {
                    my $allele2 = $ra_allele_seqs->[$j];
                    push @genotypes, "$allele1:$allele2";
                }
            }
        }
    }
 
    my $self = { reference_allele => $reference_allele,
                 coverage => $ra_coverage,
                 reference => $reference,
                 position => $position,
                 default_error => $default_error,
                 single_copy => $single_copy,
                 allele_seqs => $ra_allele_seqs,
                 genotypes => \@genotypes,
                 verbose => $verbose };

    bless $self, $class;
    $self->_check_coverage(); # make sure coverage is complete

    return $self;

} ## end new

###########################################################

=item reference_base()

  The method gets or sets the value of "reference_base", 
  which is the reference base at this position.

  Input: Optional argument sets value.
  Output: Scalar string.

=cut

###########################################################
sub reference_base {
    my $self  = shift;

    if (defined (my $new_reference_base = shift))
    {
        $self->{reference_allele} = $new_reference_base;
    }

    return $self->{reference_allele};

} ## end reference_base

###########################################################

=item reference_allele()

  The method gets or sets the value of "reference_allele", 
  which is the reference base at this position.

  Input: Optional argument sets value.
  Output: Scalar string.

=cut

###########################################################
sub reference_allele {
    my $self  = shift;

    if (defined (my $new_reference_allele = shift))
    {
        $self->{reference_allele} = $new_reference_allele;
    }

    return $self->{reference_allele};

} ## end reference_allele

###########################################################

=item coverage()

  The method returns a reference to an array of hash
  references, each of which represents an aligned short
  sequence read base at this position.  Each hash can have 
  the following fields:

      base - the base identity (required)
      error_prob - the probability that this base was
           read incorrectly (optional)
      solexa_scores - a reference to an array of four
           Solexa-style scores, for A,C,G, and T
           respectively.

  Input: Optional argument sets value.
  Output: Reference to an array of hash references.

=cut

###########################################################
sub coverage {
    my $self  = shift;
    
    if (defined (my $new_coverage = shift))
    {
        $self->{coverage} = $new_coverage;
        $self->_check_coverage();
    }

    return $self->{coverage};

} ## end coverage

###########################################################

=item coverage_string()

  The method returns a string with characters for each of
  the bases included in the coverage at this site, sorted
  by the base identities.

  Input: None.
  Output: Scalar string.

=cut

###########################################################
sub coverage_string {
    my $self  = shift;

    if (!defined ($self->{'coverage_string'}))
    {
        my $ra_coverage = $self->coverage();
        my @bases_here = map {$_->{'base'}} @{$ra_coverage};
        @bases_here = sort @bases_here;
        $self->{'coverage_string'} = join '', @bases_here;
    }

    return $self->{'coverage_string'};

} ## end coverage_string
    
###########################################################

=item add_coverage()

  The method adds the array of coverage data referenced in
  the argument to the object's current array, and returns 
  a reference to an array of hash references, each of which 
  represents an aligned short sequence read base at this 
  position (see coverage method).

  Input: Reference to an array of hash references to add to
      the current array.
  Output: Reference to an array of hash references.

=cut

###########################################################
sub add_coverage {
    my $self  = shift;
    my $ra_addl_covg = shift;
    
    push @{$self->{coverage}}, @{$ra_addl_covg};

    return $self->{coverage};

} ## end add_coverage
    
###########################################################

=item _check_coverage()

  The method checks each element in the array referenced
  by the coverage field, making sure each element has a
  "base" field, and populating the error_prob field if it
  isn't already populated.

  Input: None.
  Output: 1 if okay, croaks otherwise

=cut

###########################################################

sub _check_coverage
{
    my $self = shift;
    my $ra_coverage = $self->{coverage};

    # check that bases are populated:

    foreach my $rh_base (@{$ra_coverage})
    {
        if (!defined($rh_base->{'error_prob'}))
        {
            my $default_error = $self->{default_error};
            $rh_base->{'error_prob'} = ($rh_base->{'solexa_scores'}) ? 
                          base_error_prob($rh_base->{'solexa_scores'}) :
                          $default_error;
        }
        croak "Must define base field of each array element in -coverage passed to GTB::Var::SiteProfile constructor!\n" if (!defined ($rh_base->{'base'}));
    }

    return 1;

} ## end _check_coverage

###########################################################

=item default_error()

  The method gets or sets the value of "default_error", 
  which is the value to be used as the base_error when no
  error measure is reported.

  Input: Optional argument sets value.
  Output: Reference name (scalar string).

=cut

###########################################################
sub default_error {
    my $self  = shift;

    if (defined (my $new_default_error = shift))
    {
        $self->{default_error} = $new_default_error;
    }

    return $self->{default_error};

} ## end default_error

###########################################################

=item reference()

  The method gets or sets the value of "reference", 
  which is the name of the reference sequence the site
  resides in.

  Input: Optional argument sets value.
  Output: Reference name (scalar string).

=cut

###########################################################
sub reference {
    my $self  = shift;

    if (defined (my $new_reference = shift))
    {
        $self->{reference} = $new_reference;
    }

    return $self->{reference};

} ## end reference

###########################################################

=item position()

  The method gets or sets the value of "position", 
  which is the name of the position sequence the site
  resides in.

  Input: Optional argument sets value.
  Output: Reference name (scalar string).

=cut

###########################################################
sub position {
    my $self  = shift;

    if (defined (my $new_position = shift))
    {
        $self->{position} = $new_position;
    }

    return $self->{position};

} ## end position

###########################################################

=item filter_coverage()

  The method removes bases from the coverage array if
  requirements passed in the arguments are not fulfilled.

  Input: -min_ds_coverage will remove all bases where
           the base nucleotide is not represented at 
           least this many times on BOTH strands.
  Output: New coverage array.

=cut

###########################################################
sub filter_coverage {
    my $self  = shift;
    my %params = @_;
    my $min_ds_coverage = $params{-min_ds_coverage};

    if ($min_ds_coverage)
    {
        my $ra_coverage = $self->coverage(); 
        my $no_old_bases = @{$ra_coverage};
        my @new_coverage = ();

        foreach my $for_base (qw(A T G C))
        {
            my $rev_base = lc $for_base;
            my $no_for = grep {$_->{'base'} eq $for_base} @{$ra_coverage};
            my $no_rev = grep {$_->{'base'} eq $rev_base} @{$ra_coverage};
            next if ($no_for < $min_ds_coverage || $no_rev < $min_ds_coverage);
            push @new_coverage, 
                  grep {$_->{'base'} =~ /^($for_base|$rev_base)$/} 
                          @{$ra_coverage};
        }
        my $no_new_bases = @new_coverage;
        $self->{coverage} = \@new_coverage;
    }

    return $self->{coverage};

} ## end filter_coverage

###########################################################

=item log_likelihoods()

  The method returns a reference to a hash of log likelihood
  values, where the fields are the various genotypes (e.g.,
  'AA'), and the values are the log of the likelihood of 
  that genotype given the observed bases.

  If non-single-base alleles have been passed to the 
  constructor with the "-allele_seqs" option, the genotypes
  returned as keys to the hash will have alleles separated by
  a ":" character (since their lengths are not restricted to 1).

  Input: -prior_diff => $diff will adjust the likelihood of
      non-reference genotypes downward by $diff.
  Output: Reference to a hash of scalar numbers.

=cut

###########################################################
sub log_likelihoods {
    my $self  = shift;
    my %params = @_;
    my $prior_diff = $params{-prior_diff} || 0;
    my $align_bias = $params{-alignment_bias} || 0;

    if ($align_bias > 0.5)
    {
        croak "Cannot specify an alignment bias in log_likelihoods method that is greater than 0.5!\n";
    }

    if (!defined ($self->{log_likelihoods}))
    {
        my $ra_coverage = $self->coverage();
    
        my %log_likelihood = ();
        foreach my $genotype (@{$self->{genotypes}})
        {
            # add prior diff if necessary:
            my $ref_allele = $self->reference_allele();

            print "Evaluating genotype $genotype\n" if ($self->{verbose});
            $log_likelihood{$genotype} = 
                  (!($prior_diff) || 
                    $self->_is_homozygous_ref($genotype)) ? 0 :
                                            -1 * $prior_diff;

            foreach my $rh_base (@{$ra_coverage})
            {
                $log_likelihood{$genotype} += 
                       $self->base_log_likelihood($genotype, $rh_base, $ref_allele, $align_bias);
                print "\t$log_likelihood{$genotype}\n" if ($self->{verbose});
            }
        }

        $self->{log_likelihoods} = \%log_likelihood;
    }

    return $self->{log_likelihoods};

} ## end log_likelihoods

###########################################################

=item _is_homozygous_ref()

  This method checks a genotype against the reference allele
  and returns 1 if the genotype has only reference alleles
  present, and 0 if it has anything non-reference present.

  Input: Genotype (scalar string): colon separated if non-ATGC
      alleles and more than one copy.
  Output: 1 if homozygous ref, 0 otherwise.

=cut

###########################################################
sub _is_homozygous_ref {
    my $self = shift;
    my $genotype = shift;

    my @genotype_array = ($self->{single_copy}) ? ($genotype) :
                         ($self->{allele_seqs}) ? split /:/, $genotype :
                                                  split //, $genotype;

    my $ref_allele = $self->reference_allele();
    foreach my $allele (@genotype_array) {
        if ($allele ne $ref_allele) {
            return 0;
        }
    }

    return 1;

} ## end _is_homozygous_ref

###########################################################

=item base_log_likelihood()

  This method takes a genotype and a reference to hash of
  base information (one element of the coverage array) and
  returns the component of the log likelihood attributable
  to the base.

  Input: Genotype (scalar string), reference to a hash with
      fields 'base', and 'error_prob' or 'solexa_scores',
      optional reference base (could be undefined), and
      optional alignment bias (could be 0), reference allele,
      and align_bias value (optional)
  Output: Scalar number (log likelihood due to that base).

=cut

###########################################################
sub base_log_likelihood {
    my $self = shift;
    my $genotype = shift;
    my $rh_base = shift;
    my $ref_base = shift;
    my $align_bias = shift;

    my $this_allele = $rh_base->{'base'};
    $this_allele = uc $this_allele;
    my $error_prob = $rh_base->{'error_prob'} || $self->default_error();

    my @genotype_array = ($self->{single_copy}) ? ($genotype) :
                         ($self->{allele_seqs}) ? split /:/, $genotype :
                                                  split //, $genotype;

    if ($self->{single_copy}) { # this is easy!
        my $log_likelihood = ($genotype eq $this_allele) ? 
                             log(1 - $error_prob) : log($error_prob/3);
        return $log_likelihood;
    }

    #my ($allele1, $allele2) = split /:/, $genotype;
    my $no_alleles = @genotype_array; # for now we should only have two here
    if ($no_alleles != 2) {
        croak "Something wrong with new --single_copy option: don\'t have two alleles in default.\n";
    }
    my ($allele1, $allele2) = @genotype_array;
    if ($allele1 eq $allele2) # homozygote
    {
        my $log_likelihood = ($this_allele eq $allele1) ? 
                             log(1-$error_prob) : log($error_prob/3);
        return $log_likelihood;
    }
    else
    {
        my $ref_boost = 0; # if alignment bias, reference base has higher likelihood
        if ($ref_base && $align_bias)
        {
            $ref_boost = ($this_allele eq $ref_base) ? $align_bias : -1*$align_bias;
        }
        my $log_likelihood = (($this_allele eq $allele1) || 
                           ($this_allele eq $allele2)) ? 
                             log(0.5 + $ref_boost - $error_prob/3) : log($error_prob/3);
        return $log_likelihood;
    }

} ## end base_log_likelihood

###########################################################

=item strand_bias_score()

  This method takes a genotype as an argument, and returns
  a log score for the likelihood of seeing the observed
  strand bias if that genotype is correct.

  Input: Genotype (scalar string).  If -min_count parameter
     is passed, routine will return the minimum number of 
     instances of an allele on a single strand seen.
  Output: Scalar number (log likelihood for strand bias).

=cut

###########################################################
sub strand_bias_score {
    my $self = shift;
    my $genotype = shift ||
        croak "Must pass a genotype to strand_bias_score as an argument.\n";
    my %params = @_;
    my $min_count = $params{-min_count};

    my $ra_coverage = $self->coverage();

    my @bases = split //, $genotype;

    my %done = (); # record what's been done
    my $sb_score = 0; # additive with each base
    my $first_base = 1;

    foreach my $base (@bases)
    {
        next if ($done{$base});
        my $sum = 0;
        my @covs = grep {$_->{base} =~ /^$base$/i} @{$ra_coverage};
        my $total = @covs;
        my $forward = grep {$_->{base} eq $base} @{$ra_coverage};

        if ($min_count)
        {
            my $reverse = $total - $forward;
            my $min = ($forward < $reverse) ? $forward : $reverse;
            $sb_score = $min if ($first_base || $min <$sb_score);
            $first_base = 0;
            $done{$base} = 1;
            next;
        }

        if ($forward < $total/2)
        {
            for (my $i=0; $i<= $forward; $i++)
            {
                $sum += choose($total, $i);
            }
        }
        else
        {
            for (my $i=$forward; $i<= $total; $i++)
            {
                $sum += choose($total, $i);
            }
        }
        $sum *= 0.5**$total;
        $sum = 1e-32 if (!$sum);

        $sb_score = (log($sum) < $sb_score) ? log($sum) : $sb_score; 
        $done{$base} = 1;
    }

    return $sb_score;

} ## end strand_bias_score

###########################################################

=item minor_allele_proportion()

  This method returns the estimated proportions and 
  error estimates of all observed non-reference alleles.

  Input: Object.
  Output: Reference to an array of hash references.  Each
     hash has the following field value pairs:
     'minor_allele' => non-reference allele observed
     'minor_observed' => # of observations of minor allele
     'total_observed' => total number of observations
     'proportion' => best estimate of proportion
     'error' => best estimate of error of proportion

=cut

###########################################################
sub minor_allele_proportion {
    my $self = shift;

    my $ref_base = uc $self->reference_base();
    my $ra_coverage = $self->coverage();

    my %non_ref_counts = ();
    my %forward_counts = ();
    my %reverse_counts = ();
    my $total_bases = 0;
    my $total_forward = 0;
    my $total_reverse = 0;

    foreach my $base_cov (@{$ra_coverage})
    {
        my $base = uc $base_cov->{'base'};
        my $forward = ($base eq $base_cov->{'base'}) ? 1 : 0;
        next if ($base eq 'N');

        $total_bases++;
        $total_forward++ if ($forward);
        $total_reverse++ if (!$forward);

        if ($base ne $ref_base)
        {
            $forward_counts{$base}++ if ($forward);
            $reverse_counts{$base}++ if (!$forward);
            $non_ref_counts{$base}++;
        }
    }

    my $ra_return_props = [];
    foreach my $nr_base (keys %non_ref_counts)
    {
        my $nr_count = $non_ref_counts{$nr_base} || 0;
        my $nr_forward = $forward_counts{$nr_base} || 0;
        my $nr_reverse = $reverse_counts{$nr_base} || 0;
        my $prop = $nr_count/$total_bases;
        my $error = sqrt($prop*(1-$prop)/$total_bases);
        push @{$ra_return_props}, {'minor_allele' => $nr_base,
                                   'minor_observed' => $nr_count,
                                   'total_observed' => $total_bases,
                                   'minor_forward' => $nr_forward,
                                   'total_forward' => $total_forward,
                                   'minor_reverse' => $nr_reverse,
                                   'total_reverse' => $total_reverse,
                                   'proportion' => $prop,
                                   'error' => $error };
    }

    return $ra_return_props;
} 

###########################################################

=item one_base_ll_diff()

  This method takes the parameters "-forwards" and
  "-reverses", and reports the log likelihood score (diff
  between most likely and next most likely ) for
  the reference genotype when the base is represented
  "-forwards" times on the forward strand, and "-reverses"
  time on the reverse strand.

  Input: -forwards, -reverses, and
         -align_bias - parameter representing increased chance of seeing ref base

  Output: Scalar number (log likelihood score).

=cut

###########################################################
sub one_base_ll_diff {
    my %params = @_;
    my $forwards = $params{-forwards};
    my $reverses = $params{-reverses};
    my $align_bias = $params{-align_bias};

    my $ra_coverage = [];

    for (my $i=0; $i < $forwards; $i++)
    {
        push @{$ra_coverage}, {'base' => 'A'};
    }
    for (my $i=0; $i < $reverses; $i++)
    {
        push @{$ra_coverage}, {'base' => 'a'};
    }

    my $site_prof = GTB::Var::SiteProfile->new(
                       -reference_allele => 'A',
                       -coverage => $ra_coverage);
    my $rh_scores = $site_prof->log_likelihoods(
                       -align_bias => $align_bias);

    my @sorted_genotypes = sort {$rh_scores->{$b} <=> $rh_scores->{$a}}
                      keys %{$rh_scores};

    my $score = $rh_scores->{$sorted_genotypes[0]} - 
                         $rh_scores->{$sorted_genotypes[1]};

    return $score;

} # end one_base_ll_diff

###########################################################

sub choose
{
    my $N = shift;
    my $n = shift;

    if (($n==0) || ($n==$N))
    {
        return 1;
    }

    # do this the easier way by making $n larger:
    if ($n < $N/2)
    {
        $n = $N - $n;
    }

    my $approx = ($N - $n > 12) ? 1 : 0; # ensures 99% accuracy
    my $numerator = 1;
    if (!$approx)
    {
        for (my $i=$n+1; $i<= $N; $i++)
        {
            my $j = $i-$n;
            $numerator *= $i/$j;
        }
    }
    else
    {
        my $diff = $N - $n;
        my $log = $N*log($N) - $n*log($n) - $diff*log($diff);
        $numerator = exp($log);
        my $factor = sqrt($N/(2*3.14*$n*$diff));
        $numerator *= $factor;
    }

    return $numerator;

}

1;
__END__

=back
