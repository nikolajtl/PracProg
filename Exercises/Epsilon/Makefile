all: epsilon.exe epsilon.txt

epsilon.exe: epsilon.cc
	clang++ epsilon.cc -o epsilon.exe

epsilon.txt: epsilon.exe
	./epsilon.exe > epsilon.txt

clean:
	$(RM) epsilon.exe epsilon.txt