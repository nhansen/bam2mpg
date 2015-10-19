# $Id:$
# To create build script, run, "perl Build.PL"
use Module::Build;

my $build = Module::Build->new(
        module_name   => 'bam2mpg',
        dist_version  => '1.0.1',
        dist_author   => 'GTB Perl Programmers',
        dist_abstract => 'Script for calling variants and genotypes from single-sample BAM files',
        license       => 'perl',
        requires      => { },
        );
if (!`which samtools`) {
    warn "ERROR: samtools executable is not found.\nThis script requires "
        . "the samtools binary to be in the user's PATH\n";
}
$build->create_build_script;

