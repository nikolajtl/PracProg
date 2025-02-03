CXX = clang++
CXXFLAGS += -Wall -std=c++20 -O2
LDLIBS += -lstdc++ -lm

.PHONY: all clean run

all: math.exe Out.txt

math.exe: math.o sfuns.o
	$(CXX) $(CXXFLAGS) -o math.exe math.o sfuns.o $(LDLIBS)

math.o: math.cc sfuns.h
	$(CXX) $(CXXFLAGS) -c -o math.o math.cc

sfuns.o: sfuns.cc sfuns.h
	$(CXX) $(CXXFLAGS) -c -o sfuns.o sfuns.cc

run: math.exe
	./math.exe

Out.txt: math.exe
	./math.exe > Out.txt

clean:
	$(RM) math.exe math.o sfuns.o Out.txt