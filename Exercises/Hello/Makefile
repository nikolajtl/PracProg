Out.txt : hello.exe             # Out.txt depends on hello.exe
	./hello.exe > Out.txt        # run hello.exe, send output to Out.txt
hello.exe: hello.cc                 # hello depends on hello.cc
	clang++ -o hello.exe hello.cc  # compile-and-link hello.cc, save executable in hello.exe
clean:
	rm -f Out.txt hello.exe
