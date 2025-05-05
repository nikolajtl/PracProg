# Compiler and flags
CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

# Executable name and source files
TARGET = main.exe
SRC    = main.cc matrix.cc
OBJ    = $(SRC:.cc=.o)

OUTFILE = out.txt
HIGGSDATA = higgs_data.txt
HIGGSFIT = higgs_fit.txt

# Default rule
all: $(TARGET) run plot

# Link object files into the executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile .cc files to .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program
run: $(TARGET)
	./$(TARGET) 2>&1 | tee $(OUTFILE)

plot: $(HIGGSDATA) $(HIGGSFIT)
	gnuplot -e "\
	set terminal pdf; set output 'plot.pdf'; \
	set xlabel 'E [GeV]'; set ylabel 'σ(E)'; \
	plot 'higgs_data.txt' using 1:2:3 with yerrorbars pointtype 7 pointsize 0.5 title 'Data', \
	     'higgs_fit.txt' using 1:2 with lines linecolor rgb 'black' title 'Fit'"

# Clean up build artifacts
clean:
	rm -f $(OBJ) $(TARGET) $(OUTFILE) $(HIGGSFIT) plot.pdf

rebuild: clean all