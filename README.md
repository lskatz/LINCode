# lincodes_LK (minimal repo)

This repository contains the lincodes_LK.pl script (flat-file version) and a minimal Test Anything Protocol directory with a simple test.

What's included:
- scripts/lincodes_LK.pl — the main script (unchanged)
- test-anything-protocol/ — TAP tests (01_basic.t)
- Makefile.PL — declares CPAN prerequisites used by the script
- .github/workflows/perl.yml — GitHub Actions workflow to run tests on push / PR

CI notes
- The script depends on PDL, PDL::IO::FastRaw and File::Map which are installed in CI.
- The workflow installs prerequisites using `cpanm --installdeps .` which reads the distribution metadata / Makefile.PL.

How to run tests locally
1. Install cpanminus (cpanm) if you don't have it.
2. Install dependencies from the distribution metadata:
   - cpanm --notest --installdeps .
3. Run tests:
   - prove -l test-anything-protocol

License
- The script contains GPLv3 notices inherited from the original file. This repository is a minimal wrapper for testing and CI.
