MKDIR_P = mkdir -p
CC = gcc
LIBS = -lm
LOGLVL = 0
CFLAGS = -pg -g -Wall -I${SRC_DIR} -DLOG_LEVEL=${LOGLVL}

OBJ_DIR = ./obj
BIN_DIR = ./bin
SRC_DIR = ./src
RES_DIR = ./res
TESTS_DIR = ./tests

P_DEF_O = ${OBJ_DIR}/p-def.o
P_ARITHM_O = ${OBJ_DIR}/p-def.o ${OBJ_DIR}/p-arithm.o
P_ANALYSIS_O = ${OBJ_DIR}/p-def.o ${OBJ_DIR}/p-arithm.o ${OBJ_DIR}/p-analysis.o
P_CAUCHY_O = ${OBJ_DIR}/p-def.o ${OBJ_DIR}/p-arithm.o ${OBJ_DIR}/p-analysis.o ${OBJ_DIR}/cauchy.o

.PHONY: directories
all: directories test1 test2 test3 test4 test5 test6 test7 test8
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

test8: directories test8.o p-def.o p-arithm.o p-analysis.o cauchy.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test8.o ${P_CAUCHY_O} -o ${BIN_DIR}/test8 $(LIBS)

tar:
	tar czvf ../p-adic.tar.gz ../p-adic/Makefile ../p-adic/README ../p-adic/src ../p-adic/tests

clean:
	@rm -rf obj bin res 2&>/dev/null

