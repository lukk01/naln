#!/bin/bash
## gcc -c -g -Wall -O0 -DHAVE_PTHREAD abwape.c
## gcc -c -g -Wall -O0 -DHAVE_PTHREAD abwase.c
#gcc -c -g -Wall -O0 -DHAVE_PTHREAD abwapemain.c
gcc -c -g -Wall -O0 abwape.c
gcc -c -g -Wall -O0 abwase.c
gcc -c -g -Wall -O0 abwasemain.c
#gcc -c -g -Wall -O0 abwape.c
#gcc -c -g -Wall -O0 abwapemain.c
gcc -g -Wall -O0 -DHAVE_PTHREAD -I/Users/lukk01/tmp/naln/kustos/naln/bwape/bwa-0.6.1/ -L/Users/lukk01/tmp/naln/kustos/naln/bwape/bwa-0.6.1/ QSufSort.o bwt_gen.o utils.o bwt.o bwtio.o bwtaln.o bwtgap.o is.o bntseq.o bwtmisc.o bwtindex.o ksw.o stdaln.o simple_dp.o bwaseqio.o bwase.o kstring.o cs2nt.o bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o bwtsw2_chain.o bamlite.o fastmap.o bwtsw2_pair.o abwape.o abwase.o abwasemain.o -o abwase -lm -lz -lpthread
