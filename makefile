MAIN = test5
# other options: test3, testdr, testsw, testfa, testtools
# private options: sw2pdg, testdr2

makewhat: $(MAIN)
exe = $(MAIN)

CC = g++
CFLAGS = -g -std=c++11

.cc.o:
	$(CC) $(CFLAGS) -c $<

$(MAIN):
	$(CC) $(CFLAGS) -o $(exe) $^
	
PHONY : clean
clean : 
	rm  $(exe)

test3: test3.o classes.o ff.o EWPOZ.o xscnnlo.o B0.o C0.o D0.o li.o
test4: test4.o classes.o ff.o EWPOZ.o EWPOZ2.o xscnnlo.o B0.o C0.o D0.o li.o linex.o delrho.o
test5: test5.o classes.o ff.o EWPOZ.o EWPOZ2.o xscnnlo.o B0.o C0.o D0.o li.o linex.o delrho.o
testdr: testdr.o deltar.o delrho.o B0.o li.o linex.o
testdr2: testdr2.o deltar.o delrho.o B0.o li.o linex.o
testsw: testsw.o EWPOZ.o EWPOZ2.o delrho.o classes.o ff.o B0.o C0.o D0.o li.o linex.o
sw2pdg: sw2pdg.o EWPOZ.o EWPOZ2.o delrho.o classes.o ff.o B0.o C0.o D0.o li.o linex.o
testfa: testfa.o EWPOZ.o EWPOZ2.o delrho.o classes.o ff.o B0.o C0.o D0.o li.o linex.o
testtools: testtools.o EWPOZ.o EWPOZ2.o delrho.o classes.o ff.o B0.o C0.o D0.o li.o linex.o tools.o
