OBJS = classes.cc ff.cc FASW1.cc xscnnlo.cc B0.cc C0.cc D0.cc li.cc test1.cc
CC = g++
CFLAGS = -I
exe = output
test: $(OBJS)
	$(CC) -o $(exe) $(OBJS)


.PHONY : clean
clean : 
	rm  $(exe)
  
