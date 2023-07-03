LIB_SOURCES = $(wildcard src/*.cpp)
CFLAGS = -std=c++20 -Wall -Wextra -Werror -pedantic -O3
LIBS = -lgmp -lgmpxx

.PHONY: all
all: bin/test_pohlig_hellman bin/test_bsgs_parallel bin/test_pohlig_hellman_parallel


bin/test_pohlig_hellman: $(LIB_SOURCES) tests/test_pohlig_hellman.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_pohlig_hellman.cpp -o bin/test_pohlig_hellman $(LIBS)

bin/test_bsgs_parallel: $(LIB_SOURCES) tests/test_bsgs_parallel.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_bsgs_parallel.cpp -o bin/test_bsgs_parallel $(LIBS)

bin/test_pohlig_hellman_parallel: $(LIB_SOURCES) tests/test_pohlig_hellman_parallel.cpp
	g++ $(CFLAGS) -I./include/ -L./lib/ $(LIB_SOURCES) tests/test_pohlig_hellman_parallel.cpp -o bin/test_pohlig_hellman_parallel $(LIBS)