CC=gcc
PROG=preprocess preprocess-no-merging pileup
INDIR=src
CFLAGS=-g -Wall -O2 -Wno-unused-function

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

all:$(PROG)

preprocess:$(INDIR)/kthread.o $(INDIR)/preprocess.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

preprocess-no-merging:$(INDIR)/kthread.o $(INDIR)/preprocess-no-merging.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

pileup:$(INDIR)/kthread.o $(INDIR)/bgzf.o $(INDIR)/razf.o $(INDIR)/hts.o $(INDIR)/bedidx.o $(INDIR)/faidx.o $(INDIR)/sam.o $(INDIR)/pileup.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_MT bgzf.c -o $@

clean:
		rm -fr gmon.out $(INDIR)/*.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

preprocess.o: kvec.h kseq.h
preprocess-no-merging.o: kvec.h kseq.h
bedidx.o: ksort.h kseq.h khash.h
bgzf.o: bgzf.h
faidx.o: faidx.h khash.h razf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
pileup.o: sam.h bgzf.h hts.h faidx.h ksort.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h

