use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME         => 'lincodes_LK',
    VERSION      => '0.01',
    ABSTRACT     => 'Flat-file LINcode assignment script',
    AUTHOR       => 'Keith Jolley <keith.jolley@biology.ox.ac.uk>',
    LICENSE      => 'perl',
    PREREQ_PM    => {
        'PDL'               => 0,
        'PDL::IO::FastRaw'  => 0,
        'File::Map'         => 0,
    },
);
