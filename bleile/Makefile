# approxcolors Makefile

#
# Valgrind does not support fegetround & fesetround. With following compile option
# their use is circumvented. We also recommend to use QSopt as the LP-solver while
# debugging with valgrind, as the commercial solvers impose valgrind errors internally.
#
# CALLGRIND :
# 	valgrind --tool=callgrind ./approxcolors test/dimacs/DSJC250.5.col
#		quick inspect on https://www.speedscope.app/
#
# VALLGRIND MEMORY LEAK CHECK:
#		valgrind --leak-check=full --track-origins=yes ./approxcolors ../test/dimacs/queen10_10.col
# CFLAGS+= -DCOMPILE_FOR_VALGRIND
# CXXFLAGS += -g
#

BLEILE_DIR=.
BLEILE_HEADER=$(BLEILE_DIR)/header
BLEILE_SRC=$(BLEILE_DIR)/src
BLEILE_UTILS=$(BLEILE_DIR)/utils

# NOTE: CXXFLAGS are passed here from Makefile in parent directory
CXX=g++
CXXFLAGS += -std=c++17 -O3

all: libapproxcolors.a approxcolors

testall:
	$(foreach file, $(wildcard $(EXACTCOLOR_DIR)/test/difficult_subset/*), echo $(file);)
	$(foreach file, $(wildcard $(EXACTCOLOR_DIR)/test/difficult_subset/*), ./approxcolors $(file);)

# best results: http://cedric.cnam.fr/~porumbed/graphs/

testsingle:
	./approxcolors ../test/dimacs/flat1000_60_0.col

testdummy:
	./approxcolors ../test/dimacs/DSJC250.5.col


clean_ext:
	rm -f *.o approxcolors

clean:
	rm -f *.o *.a approxcolors


approxcolors: libapproxcolors.a $(BLEILE_DIR)/approxcolors.cpp
	$(CXX) $(CXXFLAGS) -c $(BLEILE_DIR)/approxcolors.cpp -o approxcolors.o
	$(CXX) $(CXXFLAGS) -o approxcolors approxcolors.o -L. -lapproxcolors

#  mmt.o mmt_graph.o mmt_partial_coloring.o hungarian.o

libapproxcolors.a: hungarian.o mmt_graph.o mmt_partial_coloring.o mmt.o c_connector.o
	$(AR) rcs libapproxcolors.a hungarian.o mmt_graph.o mmt_partial_coloring.o mmt.o c_connector.o

hungarian.o: $(BLEILE_UTILS)/hungarian.cpp $(BLEILE_UTILS)/hungarian.h
	$(CXX) $(CXXFLAGS) -c $(BLEILE_UTILS)/hungarian.cpp -o hungarian.o

mmt_graph.o: $(BLEILE_SRC)/mmt_graph.cpp $(BLEILE_HEADER)/mmt_graph.h
	$(CXX) $(CXXFLAGS) -c $(BLEILE_SRC)/mmt_graph.cpp -o mmt_graph.o

mmt_partial_coloring.o: hungarian.o mmt_graph.o $(BLEILE_SRC)/mmt_partial_coloring.cpp $(BLEILE_HEADER)/mmt_partial_coloring.h
	$(CXX) $(CXXFLAGS) -c $(BLEILE_SRC)/mmt_partial_coloring.cpp -o mmt_partial_coloring.o

mmt.o: mmt_graph.o mmt_partial_coloring.o $(BLEILE_SRC)/mmt.cpp $(BLEILE_HEADER)/mmt.h
	$(CXX) $(CXXFLAGS) -c $(BLEILE_SRC)/mmt.cpp -o mmt.o

c_connector.o: mmt.o mmt_graph.o $(BLEILE_SRC)/c_connector.cpp $(BLEILE_HEADER)/c_connector.h
	$(CXX) $(CXXFLAGS) -c $(BLEILE_SRC)/c_connector.cpp -o c_connector.o
