MAIN = test2

OBJS = classes.cc ff.cc EWPOs.cc xscnnlo.cc B0.cc C0.cc D0.cc li.cc $(MAIN).cc
CC = g++
CFLAGS = -I
exe = $(MAIN)
$(MAIN): $(OBJS)
	$(CC) -o $(exe) $(OBJS)


.PHONY : clean
clean : 
	rm  $(exe)
  
