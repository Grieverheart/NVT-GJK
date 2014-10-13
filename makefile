SRC=$(wildcard src/*.c)
OBJ=$(patsubst src/%.c, bin/%.o, $(SRC))
EXE=nvt-gjk

CC=gcc
CFLAGS=-Wall -O3 -std=c99 -g -msse3
LDFLAGS=-lm
dSFMTFLAGS=-msse2 -DDSFMT_MEXP=19937 -DHAVE_SSE2
RM=rm

bin/%.o: src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

bin/dSFMT.o: src/dSFMT.c
	$(CC) $(CFLAGS) $(dSFMTFLAGS) -o $@ -c $<
	
bin/boop.o: src/boop.c include/boop.h
	$(CC) $(CFLAGS) -o $@ -c $<

bin/quaternion.o: src/quaternion.c include/quaternion.h
	$(CC) $(CFLAGS) -o $@ -c $<

bin/gjk.o: src/gjk.c include/gjk.h
	$(CC) $(CFLAGS) -o $@ -c $<

bin/nvt_gjk.o: src/nvt_gjk.c
	$(CC) $(CFLAGS) -o $@ -c $<
	
.PHONY: all
all: $(EXE)
	@echo Done

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)
	
.PHONY: clean
clean:
	-$(RM) $(OBJ)
	@echo Clean Done!
