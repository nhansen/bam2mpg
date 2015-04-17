# $Id$
# t/02_call.t - check that shimmer successfully calls variants

use strict;
use Test::More;
use Module::Build;

my $bam2mpg = 'scripts/bam2mpg';

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => 3;

my $out;
my $testref = 't/testref.fa';
my $testbam = 't/testbam.bam';
my $testbedfile = 't/testbed.bed';

system("perl -w -I lib $bam2mpg $testref $testbam -r chr4:488000-495000 > t/calltest1.out 2>&1");
$out = `awk '\$1=="MPG_SNV" && \$3==491041 {print \$5}' t/calltest1.out`;
like $out, qr/CG/, "$bam2mpg variant";
$out = `grep '489247:489258' t/calltest1.out`;
like $out, qr/TTTTTTTTTT/, "$bam2mpg chr4indel";
system("perl -w -I lib $bam2mpg $testref $testbam -r chr17:690000-691000 > t/calltest2.out 2>&1");
$out = `grep '690303:690306' t/calltest2.out`;
like $out, qr/\s\*:AA\s/, "$bam2mpg chr17indel";
