# Compiler and flags
CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

# Targets and source files
TARGET   = LinEq.exe
SRC      = Linear_Equations.cc

# Timing executable
TIMING_TARGET = Timing.exe
TIMING_SRC    = Timing.cc

# Linspace-like range:
START    = 100      # Start of the range
END      = 2000     # End of the range
STEP     = 100       # Step size

# Default target: build both executables
all: $(TARGET) $(TIMING_TARGET) run timing plot

# Compile LinEq.exe
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Compile Timing.exe
$(TIMING_TARGET): $(TIMING_SRC)
	$(CXX) $(CXXFLAGS) -o $(TIMING_TARGET) $(TIMING_SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

# Run timing tests with linspace-like range and save to timing.txt
timing: $(TIMING_TARGET)
	@echo "Timing the QR-factorization for different matrix sizes:"
	@echo "# N Time(s)" > timing.txt
	@for n in $$(seq $(START) $(STEP) $(END)); do \
		echo "Running for N = $$n..."; \
		/usr/bin/time -f "$$n %e" ./$(TIMING_TARGET) $$n 2>> timing.txt; \
	done
	@echo "Results saved in timing.txt"

plot: timing.txt
	gnuplot -e "\
	set terminal pdf; set output 'timing_plot.pdf'; \
	set xlabel 'N [-]'; set ylabel 'Time [s]'; \
	set xrange [0 : $(END)*1.1]; \
	set yrange [-1 : 25]; \
	f(x) = a * x**3; \
	a = 1; \
	fit f(x) 'timing.txt' using 1:2 via a; \
	print 'Fitted coefficient a = ', a; \
	set label sprintf('Fitted coefficient a = %.6e', a) at graph 0.05, 0.9; \
	set key bottom right; \
	plot 'timing.txt' with points title 'Data Points', \
	     f(x) title 'a * N^3';"



# Clean build artifacts
clean:
	rm -f $(TARGET) $(TIMING_TARGET) timing.txt timing_plot.pdf fit.log