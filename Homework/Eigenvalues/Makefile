# Meste kopieret fra Linear_Equations Makefile'en

CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET1   = eigenval.exe
SRC1      = eigenval.cc

TARGET2	  = hydrogen.exe
SRC2	  = hydrogen.cc

TARGET3	  = convergence.exe
SRC3	  = convergence.cc

all: $(TARGET1) $(TARGET2) $(TARGET3) run plotfuncs plotconv_delta plotconv_max

$(TARGET1): $(SRC1)
	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(SRC1) $(LDFLAGS)

$(TARGET2): $(SRC2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SRC2) $(LDFLAGS)

$(TARGET3): $(SRC3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SRC3) $(LDFLAGS)

run: $(TARGET1) $(TARGET3) $(TARGET3)
	./$(TARGET1) 8
	./$(TARGET2)
	./$(TARGET3)

plotfuncs: wavefuncs.txt
	gnuplot -e "\
	set terminal pdfcairo enhanced; set output 'plot_funcs.pdf'; \
	set xlabel 'r [a_0]'; set ylabel '𝛙 [1/a_0]'; \
	set xrange [-0.1 : 20]; set yrange [-0.4 : 0.9]; \
	plot for [i=2 : *] 'wavefuncs.txt' using 1:i with points title columnheader(i) "

plotconv_delta: conv_delta.txt
	gnuplot -e "\
	set terminal pdfcairo enhanced; set output 'plot_delta.pdf'; \
	set xlabel 'Δr [a_0]'; set ylabel 'ε_0'; \
	plot 'conv_delta.txt' using 1:2 with points title 'Ground state energies w.r.t. Δr' "

plotconv_max: conv_max.txt
	gnuplot -e "\
	set terminal pdfcairo enhanced; set output 'plot_max.pdf'; \
	set xlabel 'r max [a_0]'; set ylabel 'ε_0'; \
	plot 'conv_max.txt' using 1:2 with points title 'Ground state energies w.r.t. r max' "

clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) plot_funcs.pdf plot_delta.pdf plot_max.pdf wavefuncs.txt conv_delta.txt conv_max.txt