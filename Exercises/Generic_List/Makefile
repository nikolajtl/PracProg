CXX = clang++
CXXFLAGS += -Wall -std=c++20 -O2

.PHONY: all clean run

# Default target to compile the program
all: genlist.exe run

# Rule to create the executable from main.o
genlist.exe: main.o
	$(CXX) $(CXXFLAGS) -o genlist.exe main.o

# Rule to compile main.cc into main.o
main.o: main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

# Rule to clean up generated files
clean:
	$(RM) genlist.exe main.o output.txt

# Run the program with input redirection and output redirection
run: genlist.exe
	./genlist.exe < inputfile.txt > output.txt
