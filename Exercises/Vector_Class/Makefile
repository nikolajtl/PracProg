all: vec.exe vector_class.so vec.txt

vec.exe: main.cc vec_class.cc vec_class.h
	clang++ main.cc vec_class.cc -o vec.exe

vector_class.so: vec_class.cc
	clang++ -shared -fPIC vec_class.cc -o vector_class.so

vec.txt: vec.exe
	./vec.exe > vec.txt

clean:
	rm -f vec.exe vec.txt vector_class.so
