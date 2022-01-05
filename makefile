OBJS = classes.cc ff.cc B0.cc C0.cc D0.cc li.cc main.cc
CC = g++
CFLAGS = -I
box: $(OBJS)
	$(CC) -o ogtest $(OBJS)


.PHONY : clean
clean : 
	rm  *.o 

  