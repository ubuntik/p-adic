CXX = gcc
CFLAGS = -g -Wall -std=c99
LIBS = -lm
OBJS = p-adic.o
SRCS = p-adic.c test1.c test2.c test3.c
HDRS = p-adic.h

all: test1 test2 test3 test4 test5 test6

test1: test1.o p-adic.o
	$(CXX) $(CFLAGS) test1.o p-adic.o -o test1 $(LIBS)

test2: test2.o p-adic.o
	$(CXX) $(CFLAGS) test2.o p-adic.o -o test2 $(LIBS)

test3: test3.o p-adic.o
	$(CXX) $(CFLAGS) test3.o p-adic.o -o test3 $(LIBS)

test4: test4.o p-adic.o
	$(CXX) $(CFLAGS) test4.o p-adic.o -o test4 $(LIBS)

test5: test5.o p-adic.o
	$(CXX) $(CFLAGS) test5.o p-adic.o -o test5 $(LIBS)

test6: test6.o p-adic.o
	$(CXX) $(CFLAGS) test6.o p-adic.o -o test6 $(LIBS)

test1.o: p-adic.h test1.c
	$(CXX) $(CFLAGS) -c test1.c -o test1.o

test2.o: p-adic.h test2.c
	$(CXX) $(CFLAGS) -c test2.c -o test2.o

test3.o: p-adic.h test3.c
	$(CXX) $(CFLAGS) -c test3.c -o test3.o

test4.o: p-adic.h test4.c
	$(CXX) $(CFLAGS) -c test4.c -o test4.o

test5.o: p-adic.h test5.c
	$(CXX) $(CFLAGS) -c test5.c -o test5.o

test6.o: p-adic.h test6.c
	$(CXX) $(CFLAGS) -c test6.c -o test6.o

p-adic.o: p-adic.h p-adic.c
	$(CXX) $(CFLAGS) -c p-adic.c -o p-adic.o

clean:
	rm *.o test1 test2 test3 test4 test5 test6

