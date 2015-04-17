# $Id:$
# To create build script, run, "perl Build.PL"
use Module::Build;

my $class = Module::Build->subclass(
	class => 'My::Build',
	code => q{
		sub ACTION_code {
			my $self = shift;
                        $self->SUPER::ACTION_code();
			$self->do_system(qw(make));
		}
		} );

my $build = $class->new(
        module_name   => 'bam2mpg',
        dist_version  => 1.0,
        dist_author   => 'GTB Perl Programmers',
        dist_abstract => 'Script for calling variants and genotypes from single-sample BAM files',
        license       => 'perl',
        scripts       => [qw(bam2mpg)],
        requires      => { },
        );
if (!`which samtools`) {
    warn "ERROR: samtools executable is not found.\nThis script requires "
        . "the samtools binary to be in the user's PATH\n";
}
$build->create_build_script;

