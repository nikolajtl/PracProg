# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O2 -pthread

# Source and executable names
SRC = main.cc
EXE = main.exe

# Time command and output file
TIME = time --portability --append --output $@

# Define the number of terms for the calculation
N = 1e9

# Default target
all: $(EXE) Out.txt

# Rule to build the executable
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(SRC)

# Target to run the program with different thread counts and log output to Out.txt
Out.txt: $(EXE) Makefile
	>$@
	$(TIME) ./$(EXE) -terms $(N) -threads 1 >> $@
	$(TIME) ./$(EXE) -terms $(N) -threads 2 >> $@
	$(TIME) ./$(EXE) -terms $(N) -threads 3 >> $@
	$(TIME) ./$(EXE) -terms $(N) -threads 4 >> $@

# Clean up generated files
clean:
	rm -f $(EXE) Out.txt

# Run the executable (you can add parameters here if needed)
run: $(EXE)
	./$(EXE)

.PHONY: all clean run
