CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

SRC = main.cc matrix.cc
HDR = matrix.h
OUT = main.exe

HYDROGEN_SRC = hydrogen.cc matrix.cc
HYDROGEN_OUT = hydrogen.exe

JACOBI_SRC = jacobi.cc matrix.cc
JACOBI_OUT = jacobi.exe

all: $(OUT) $(HYDROGEN_OUT) $(JACOBI_OUT) run plot

# Linking for main.cc
$(OUT): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(SRC)

# Linking for hydrogen.cc
$(HYDROGEN_OUT): $(HYDROGEN_SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $(HYDROGEN_OUT) $(HYDROGEN_SRC)

# Linking for jacobi.cc
$(JACOBI_OUT): $(JACOBI_SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $(JACOBI_OUT) $(JACOBI_SRC)

run: $(OUT) $(HYDROGEN_OUT)
	./$(OUT) > out_main.txt
	./$(HYDROGEN_OUT) > out_hydrogen.txt
	./$(JACOBI_OUT) > out_jacobi.txt

plot:
	gnuplot -persist -e "\
		set terminal pdf; set output 'convplot.pdf'; \
		set title 'Convergence w.r.t. n, N=200'; \
		set xlabel 'n'; \
		set ylabel 'Ground state energy'; \
		set grid; \
		plot 'convergence.txt' using 1:2 with linespoints title ' ';"

clean:
	rm -f $(OUT) $(HYDROGEN_OUT) $(JACOBI_OUT) out_main.txt out_hydrogen.txt out_jacobi.txt convergence.txt convplot.pdf

rebuild: clean all