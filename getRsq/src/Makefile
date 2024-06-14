
CC=gcc
CFLAGS := ${CFLAGS} -O3 #-std=gnu99 -g -O2
LDFLAGS := ${LDFLAGS} -lhts -lgsl -lgslcblas -llapack -lm -lz #-static -fno-lto #-fcommon #-static

PROGRAM = getRsq
PROGFLAG = -DGETRSQ
OBJS = loadVCF.o getLogBF.o sort.o util.o getVCFHeader.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LDFLAGS) $(PROGFLAG)

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)

