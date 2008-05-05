CC		= g++
CLINKER		= g++
CFLAGS		=
LIBS		= -lm
OPTF		= -O

#Final output name of the execution file
TARGETS 	= gtest

all:      $(TARGETS)

OBJS		= gtest.o gport.o

$(TARGETS): $(OBJS) 
	$(CC) $(CFLAGS) $(OPTF) -o $(TARGETS) $(OBJS) $(LIBS)


gtest.o: gtest.cpp
	$(CC) $(CFLAGS) $(OPTF) -c gtest.cpp 

gport.o: gport.h gport.cpp
	$(CC) $(CFLAGS) $(OPTF) -c gport.cpp  

clean:
	rm -f *.o *.bak core
