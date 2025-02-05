all: vec.exe vec.txt

vec.exe: vec_main.cc vec_class.cc vec_class.h
	clang++ vec_main.cc vec_class.cc -o vec.exe

vec.txt: vec.exe
	./vec.exe > vec.txt

clean:
	$(RM) vec.exe vec.txt