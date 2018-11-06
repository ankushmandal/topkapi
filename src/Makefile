CC=mpiCC
STDFLAG=-std=gnu++11
OFLAGS=-O3
DFLAGS=-D__USE_LARGEFILE64 -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
OMPFLAG=-fopenmp
CFLAGS=${STDFLAG} ${OFLAGS} ${DFLAGS} ${OMPFLAG}

DEPS=MurmurHash3.cpp HashFunction.cpp LossyCountMinSketch.cpp
MAIN=topkapi.cpp
SOURCES=$(DEPS) $(MAIN)
HEADERS=$(DEPS:.cpp=.h)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=$(basename $(MAIN))

all: $(SOURCES) $(HEADERS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OMPFLAG) $(OBJECTS) -o $@

$(DEPS:.cpp=.o): %.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $< -o $@

$(MAIN:.cpp=.o): %.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJECTS)

.PHONY: cleaner

cleaner: clean
	rm -f $(EXECUTABLE)
