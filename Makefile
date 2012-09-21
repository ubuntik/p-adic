MKDIR_P = mkdir -p
CC = gcc
CFLAGS = -g -Wall -std=c99
LIBS = -lm
OBJ_DIR = ./obj
BIN_DIR = ./bin
SRC_DIR = ./src
TESTS_DIR = ./tests

.PHONY: directories
all: directories test1 test2 test3 test4 test5 test6
directories: ${OBJ_DIR} ${BIN_DIR} ${SRC_DIR} ${TESTS_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

%.o: ${SRC_DIR}/%.c ${SRC_DIR}/p-adic.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o ${OBJ_DIR}/$@

%.o: ${TESTS_DIR}/%.c ${SRC_DIR}/p-adic.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o ${OBJ_DIR}/$@

test1: directories test1.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test1.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test1 $(LIBS)

test2: directories test2.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test2.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test2 $(LIBS)

test3: directories test3.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test3.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test3 $(LIBS)

test4: directories test4.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test4.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test4 $(LIBS)

test5: directories test5.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test5.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test5 $(LIBS)

test6: directories test6.o p-adic.o
	$(CC) $(CFLAGS) ${OBJ_DIR}/test6.o ${OBJ_DIR}/p-adic.o -o ${BIN_DIR}/test6 $(LIBS)

tar:
	tar czvf ../p-adic.tar.gz ../p-adic/Makefile ../p-adic/README ../p-adic/src ../p-adic/tests

clean:
	@rm -rf obj bin 2&>/dev/null

