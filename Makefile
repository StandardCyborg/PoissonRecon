PR_TARGET=PoissonRecon
PR_SOURCE=CmdLineParser.cpp MemoryFileSystem.cpp Factor.cpp Geometry.cpp MarchingCubes.cpp PlyFile.cpp PoissonRecon.cpp

CFLAGS += -fopenmp -Wno-deprecated -Wno-write-strings -std=c++11
LFLAGS += -lgomp -lstdc++

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3 

SRC = Src/
BIN = Bin/Linux/
INCLUDE = /usr/include/

CC=gcc
CXX=g++
MD=mkdir

PR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PR_SOURCE))))


all: CFLAGS += $(CFLAGS_DEBUG)
all: LFLAGS += $(LFLAGS_DEBUG)
all: $(BIN)
all: $(BIN)$(PR_TARGET)

release: CFLAGS += $(CFLAGS_RELEASE)
release: LFLAGS += $(LFLAGS_RELEASE)
release: $(BIN)
release: $(BIN)$(PR_TARGET)

clean:
	rm -f $(BIN)$(PR_TARGET)
	rm -f $(PR_OBJECTS)

$(BIN):
	$(MD) -p $(BIN)

$(BIN)$(PR_TARGET): $(PR_OBJECTS)
	$(CXX) -o $@ $(PR_OBJECTS) $(LFLAGS)

$(BIN)%.o: $(SRC)%.c
	mkdir -p $(BIN)
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	mkdir -p $(BIN)
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

