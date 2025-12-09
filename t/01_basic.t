#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use File::Spec;
use File::Basename qw/dirname/;

my $pwd = dirname($0);
my $script = "$pwd/../scripts/lincodes_LK.pl";
ok( -e $script, 'script exists' );

