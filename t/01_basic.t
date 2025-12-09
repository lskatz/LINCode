#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 1;
use File::Spec;
use File::Basename qw/dirname/;

my $pwd = dirname($0);
my $script = "$pwd/../scripts/lincodes_LK.pl";

# Make a dummy test for now but we will need to build this up
is(1, 1, "Does 1==1? Probably yes");

