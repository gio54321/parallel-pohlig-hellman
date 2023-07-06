LIB_SOURCES = $(wildcard src/*.cpp)
CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3
LIBS = -lgmp -lgmpxx

FF_INCLUDE = ../fastflow/

.PHONY: all
all: bin/test_bsgs_parallel bin/test_pohlig_hellman_parallel bin/bsgs_ff

bin/test_bsgs_parallel: $(LIB_SOURCES) tests/test_bsgs_parallel.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_bsgs_parallel.cpp -o bin/test_bsgs_parallel $(LIBS)

bin/test_pohlig_hellman_parallel: $(LIB_SOURCES) tests/test_pohlig_hellman_parallel.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_pohlig_hellman_parallel.cpp -o bin/test_pohlig_hellman_parallel $(LIBS)

bin/pohlig_hellman_ff: $(LIB_SOURCES) ff/pohlig_hellman_ff.cpp
	g++ $(CFLAGS) -I./include/ -I${FF_INCLUDE} -L./lib/ ff/pohlig_hellman_ff.cpp -o bin/pohlig_hellman_ff $(LIBS) -DBLOCKING_MODE