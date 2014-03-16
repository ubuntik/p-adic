MKDIR_P = mkdir -p
CC = gcc
CFLAGS = -pg -g -Wall -I${SRC_DIR}
LIBS = -lm

OBJ_DIR = ./obj
BIN_DIR = ./bin
SRC_DIR = ./src
RES_DIR = ./res
TESTS_DIR = ./tests

P_DEF_O = ${OBJ_DIR}/p-def.o
P_ARITHM_O = ${OBJ_DIR}/p-def.o ${OBJ_DIR}/p-arithm.o
P_ANALYSIS_O = ${OBJ_DIR}/p-def.o ${OBJ_DIR}/p-arithm.o ${OBJ_DIR}/p-analysis.o

##DEBUG ?= 1
##ifeq (DEBUG, 1)
##	CFLAGS =-g3 -gdwarf2 -DDEBUG
##else
##	CFLAGS=-DNDEBUG
##endif


.PHONY: directories
all: directories test1 test2 test3 test4 test5 test6 test7
directories: ${OBJ_DIR} ${BIN_DIR} ${SRC_DIR} ${TESTS_DIR} ${RES_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

${RES_DIR}:
	${MKDIR_P} ${RES_DIR}

# To obtain object files
%.o: ${SRC_DIR}/%.c
	$(CC) -c $(CFLAGS) $< -o ${OBJ_DIR}/$@

%.o: ${TESTS_DIR}/%.c
	$(CC) -c $(CFLAGS) $< -o ${OBJ_DIR}/$@

#cauchy.o: ${SRC_DIR}/cauchy_problem.c ${SRC_DIR}/cauchy_problem.h ${SRC_DIR}/p-analysis.h
#	$(CC) -c $(CFLAGS) ${SRC_DIR}/cauchy_problem.c -o ${OBJ_DIR}/cauchy.o

test1: directories test1.o p-def.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test1.o ${P_DEF_O} -o ${BIN_DIR}/test1 $(LIBS)

test2: directories test2.o p-def.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test2.o ${P_DEF_O} -o ${BIN_DIR}/test2 $(LIBS)

test3: directories test3.o p-def.o p-arithm.o p-analysis.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test3.o ${P_ANALYSIS_O} -o ${BIN_DIR}/test3 $(LIBS)

test4: directories test4.o p-def.o p-arithm.o p-analysis.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test4.o ${P_ANALYSIS_O} -o ${BIN_DIR}/test4 $(LIBS)

test5: directories test5.o p-def.o p-arithm.o p-analysis.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test5.o ${P_ANALYSIS_O} -o ${BIN_DIR}/test5 $(LIBS)

test6: directories test6.o p-def.o p-arithm.o p-analysis.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test6.o ${P_ANALYSIS_O} -o ${BIN_DIR}/test6 $(LIBS)

test7: directories test7.o p-def.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test7.o ${P_DEF_O} -o ${BIN_DIR}/test7 $(LIBS)

#test9: directories test9.o cauchy.o p-adic.o
#	$(CC) $(CFLAGS) ${OBJ_DIR}/test9.o ${OBJ_DIR}/p-adic.o ${OBJ_DIR}/cauchy.o -o ${BIN_DIR}/test9 $(LIBS)

tar:
	tar czvf ../p-adic.tar.gz ../p-adic/Makefile ../p-adic/README ../p-adic/src ../p-adic/tests

clean:
	@rm -rf obj bin res 2&>/dev/null

