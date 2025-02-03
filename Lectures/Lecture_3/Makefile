CXXFLAGS += -Wall -Wextra -Werror -std=c++23 -O2
LDLIBS += -lstdc++ -lm
.PHONEY: all
all: Out.txt Out++.txt

Out.txt: main.exe
	mono $< > $@ # mono main.exe > Out.txt

Out++.txt: main
	./main > Out++.txt

main:  main.o
	$(CXX) $(LDFALGS) $(LDLIBS) -o main main.o

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o main.cc

main.exe: main.cs sfuns.dll
	mcs -target:exe -reference:sfuns.dll -out:main.exe main.cs

sfuns.dll: sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs

.PHONEY: clean
clean:
	$(RM) sfuns.dll main.exe Out.txt

test:
	@echo $(RM)
	@echo $(CXX)
