MAIN = testtools
# other options: testmatel, testdeltar, testtools

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

testmatel: testmatel.o classes.o ff.o EWPOZ.o EWPOZ2.o xscnnlo.o B0.o C0.o D0.o li.o linex.o delrho.o
testdeltar: testdeltar.o deltar.o delrho.o classes.o B0.o li.o linex.o
testtools: testtools.o EWPOZ.o EWPOZ2.o delrho.o classes.o ff.o B0.o C0.o D0.o li.o linex.o tools.o
