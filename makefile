OBJS = classes.cc ff.cc FASW1.cc xscnnlo.cc B0.cc C0.cc D0.cc li.cc test1.cc
CC = g++
CFLAGS = -I
box: $(OBJS)
	$(CC) -o ogtest $(OBJS)


.PHONY : clean
clean : 
	rm  *.o 

  
