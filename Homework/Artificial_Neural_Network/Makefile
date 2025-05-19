CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET = main.exe
OBJS = main.o minimizer.o matrix.o
OUTFILE = out.txt
FITFILE = values.txt
PLOTS = fitplot.pdf intplot.pdf derivplot.pdf allplot.pdf diffeq_plot.pdf

DIFFEQ_TARGET = diffeq.exe
DIFFEQ_OBJS = diffeq.o minimizer.o matrix.o
DIFFEQ_FITFILE = diffeq_values.txt

all: $(TARGET) diffeq run plots

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

diffeq: $(DIFFEQ_TARGET)

$(DIFFEQ_TARGET): $(DIFFEQ_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

diffeq.o: diffeq.cc minimizer.h matrix.h
	$(CXX) $(CXXFLAGS) -c diffeq.cc

main.o: main.cc minimizer.h matrix.h
	$(CXX) $(CXXFLAGS) -c main.cc

minimizer.o: minimizer.cc minimizer.h matrix.h
	$(CXX) $(CXXFLAGS) -c minimizer.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

run: $(TARGET) $(DIFFEQ_TARGET)
	./$(TARGET)
	./$(DIFFEQ_TARGET)

plots: plot_fit plot_int plot_deriv plot_all plot_diffeq

plot_fit: $(FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'fitplot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	plot 'values.txt' using 1:2 with lines title 'g(x)', \
	     'values.txt' using 1:3 with lines title 'Trained response'"

plot_int: $(FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'intplot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	plot 'values.txt' using 1:2 with lines title 'g(x)', \
	     'values.txt' using 1:4 with lines title 'Trained integral'"

plot_deriv: $(FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'derivplot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	plot 'values.txt' using 1:2 with lines title 'g(x)', \
	     'values.txt' using 1:5 with lines title 'Trained derivative', \
		 'values.txt' using 1:6 with lines title 'Trained 2nd derivative'"

plot_all: $(FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'allplot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	plot 'values.txt' using 1:2 with lines title 'g(x)', \
		 'values.txt' using 1:3 with lines title 'Trained response', \
		 'values.txt' using 1:4 with lines title 'Trained integral', \
	     'values.txt' using 1:5 with lines title 'Trained derivative', \
		 'values.txt' using 1:6 with lines title 'Trained 2nd derivative'"

plot_diffeq: $(DIFFEQ_FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'diffeq_plot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	plot 'diffeq_values.txt' using 1:2 with lines title 'Trained response', \
	     'diffeq_values.txt' using 1:3 with lines title 'Analytic solution'"


clean:
	rm -f $(OBJS) $(TARGET) $(OUTFILE) $(FITFILE) $(PLOTS) $(DIFFEQ_TARGET) $(DIFFEQ_OBJS) $(DIFFEQ_FITFILE)

rebuild: clean all

.PHONY: diffeq