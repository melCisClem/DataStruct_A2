CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic

OBJS = simplify.o polygon.o geometry.o simplifier.o

simplify: $(OBJS)
	$(CXX) $(CXXFLAGS) -o simplify $(OBJS)

simplify.o: simplify.cpp polygon.h simplifier.h
	$(CXX) $(CXXFLAGS) -c simplify.cpp

polygon.o: polygon.cpp polygon.h geometry.h
	$(CXX) $(CXXFLAGS) -c polygon.cpp

geometry.o: geometry.cpp geometry.h
	$(CXX) $(CXXFLAGS) -c geometry.cpp

simplifier.o: simplifier.cpp simplifier.h polygon.h geometry.h
	$(CXX) $(CXXFLAGS) -c simplifier.cpp

clean:
	rm -f *.o simplify