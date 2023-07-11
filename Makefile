LIB_SOURCES = $(wildcard src/*.cpp)
CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3
LIBS = -lgmp -lgmpxx

FF_INCLUDE = ../fastflow/

.PHONY: all
all: bin/test_bsgs bin/test_pohlig_hellman bin/pohlig_hellman_ff

.PHONY: clean
clean:
	rm -f bin/*

bin/test_bsgs: $(LIB_SOURCES) tests/test_bsgs.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_bsgs.cpp -o bin/test_bsgs $(LIBS)

bin/test_pohlig_hellman: $(LIB_SOURCES) tests/test_pohlig_hellman.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_pohlig_hellman.cpp -o bin/test_pohlig_hellman $(LIBS)

bin/pohlig_hellman_ff: ff/pohlig_hellman_ff.cpp
	g++ $(CFLAGS) -I./include/ -I${FF_INCLUDE} -L./lib/ ff/pohlig_hellman_ff.cpp -o bin/pohlig_hellman_ff $(LIBS) -DBLOCKING_MODE