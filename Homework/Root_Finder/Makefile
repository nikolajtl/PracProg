CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

# Source files
COMMON_SRC = matrix.cc ode.cc
MAIN_SRC   = main.cc
CONV_SRC   = conv.cc

# Object files
COMMON_OBJ = $(COMMON_SRC:.cc=.o)
MAIN_OBJ   = $(MAIN_SRC:.cc=.o) $(COMMON_OBJ)
CONV_OBJ   = $(CONV_SRC:.cc=.o) $(COMMON_OBJ)

# Executables
TARGETS = main.exe conv.exe

# Output and plot
OUTFILE = out.txt
DATAFILES = shooting.txt exact.txt rmax.txt rmin.txt acceps.txt
PLOTS = plot_main.pdf plot_rmax.pdf plot_rmin.pdf plot_acceps.pdf
E := $(shell [ -f shooting.txt ] && head -n1 shooting.txt || echo "-1.0")

all: $(TARGETS) run plots

main.exe: $(MAIN_OBJ)
	$(CXX) $(MAIN_OBJ) -o $@ $(LDFLAGS)

conv.exe: $(CONV_OBJ)
	$(CXX) $(CONV_OBJ) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: main.exe conv.exe
	./main.exe 2>&1 | tee $(OUTFILE)
	./conv.exe

plots : plot_main plot_rmax plot_rmin plot_acceps

plot_main: $(DATAFILES) $(OUTFILE)
	gnuplot -e "\
	set terminal pdf; set output 'plot_main.pdf'; \
	set xlabel 'r'; set ylabel 'f(r)'; \
	set xrange [-0.1 : 8]; set yrange [-0.1 : 0.5]; \
	plot 'shooting.txt' using 1:2 with lines linestyle 2 title 'Shooting (E = $(E))', \
	     'exact.txt' using 1:2 with lines linecolor rgb 'black' title 'Exact (E = -0.5)'"

plot_rmax: $(DATAFILES)
	echo "\
	set terminal pdf; \
	set output 'plot_rmax.pdf'; \
	set xlabel 'r'; set ylabel 'f(r)'; \
	set xrange [-0.1 : 7.5]; set yrange [-0.1 : 0.5]; \
	plot for [i=0:5] 'rmax.txt' index i using 1:2 with lines dashtype (20,5) linetype i+1 title sprintf('Rmax = %d', i+2), \
	     'exact.txt' using 1:2 with lines linecolor rgb 'black' title 'Analytic'" > plot_rmax.gp && \
	gnuplot plot_rmax.gp && rm plot_rmax.gp

plot_rmin: $(DATAFILES)
	echo "\
	set terminal pdf; \
	set output 'plot_rmin.pdf'; \
	set xlabel 'r'; set ylabel 'f(r)'; \
	set xrange [-0.1 : 8.0]; set yrange [-0.1 : 0.5]; \
	plot for [i=0:4] 'rmin.txt' index i using 1:2 with lines dashtype (20,5) linetype i+1 title sprintf('Rmin = 0.%d', i+1), \
	     'exact.txt' using 1:2 with lines linecolor rgb 'black' title 'Analytic'" > plot_rmax.gp && \
	gnuplot plot_rmax.gp && rm plot_rmax.gp

plot_acceps: $(DATAFILES)
	echo "\
	set terminal pdf; \
	set output 'plot_acceps.pdf'; \
	set xlabel 'r'; set ylabel 'f(r)'; \
	set xrange [-0.1 : 8.0]; set yrange [-0.1 : 0.5]; \
	plot for [i=0:5] 'acceps.txt' index i using 1:2 with lines dashtype (20,5) linetype i+1 title sprintf('acc=eps = 1e-%d', i), \
	     'exact.txt' using 1:2 with lines linecolor rgb 'black' title 'Analytic'" > plot_rmax.gp && \
	gnuplot plot_rmax.gp && rm plot_rmax.gp

clean:
	rm -f *.o $(TARGETS) $(OUTFILE) $(DATAFILES) $(PLOTS)

rebuild: clean all
