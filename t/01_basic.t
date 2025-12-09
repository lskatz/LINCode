#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use File::Spec;

my $script = File::Spec->catfile('scripts', 'lincodes_LK.pl');
ok( -e $script, 'script exists' );

# Syntax check
my $c = `perl -c $script 2>&1`;
like( $c, qr/syntax OK/, 'perl -c reports syntax OK' );

# --help should print usage (exit 0)
my $help = `perl $script --help 2>&1`;
ok( $help =~ /Usage: lincodes_LK.pl/, '--help prints usage' );
