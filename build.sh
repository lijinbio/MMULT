#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
set -e
condaroot=$(conda info --base)
source "$condaroot/etc/profile.d/conda.sh"
conda init
conda activate
autoreconf -i
./configure CPPFLAGS=-I"$condaroot/include" LDFLAGS=-L"$condaroot/lib -Wl,-rpath,$condaroot/lib"
make
