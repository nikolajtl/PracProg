CXX = clang++
CXXFLAGS += -Wall -std=c++20 -O2

.PHONY: all clean

all: genlist.exe

genlist.o: genlist.cc genlist.h
	$(CXX) $(CXXFLAGS) -c -o genlist.o genlist.cc

genlist.exe: genlist.o
	$(CXX) $(CXXFLAGS) -o genlist.exe genlist.o

clean:
	$(RM) genlist.exe genlist.o
