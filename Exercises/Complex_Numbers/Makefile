CXX = clang++
CXXFLAGS = -std=c++20
RM = rm -f

# .PHONY to specify that these are not actual files
.PHONY: all run clean

# Default target
all: main.exe Out.txt

# Rule to build the executable
main.exe: main.cc
	$(CXX) $(CXXFLAGS) -o $@ $<

# Rule to generate Out.txt
Out.txt: main.exe
	./main.exe > Out.txt

# Rule to run the program (run the executable)
run: main.exe
	./main.exe

# Rule to clean generated files
clean:
	$(RM) main.exe Out.txt
