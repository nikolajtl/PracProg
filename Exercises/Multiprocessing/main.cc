#include<iostream>
#include<thread>
#include<string>
#include<vector>

struct datum {long start,stop; double sum;};

void harm(datum& d){
	d.sum=0;
	for (long i=d.start+1;i<=d.stop;i++)d.sum+=1.0/i;
	}

int main (int argc, char** argv) {
	long nthreads=1, nterms=(long)1e8;
	for (int i=0;i<argc;i++) {
		std::string arg = argv[i];
		if(arg=="-threads" && i+1<argc) nthreads=std::stoi(argv[i+1]);
		if(arg=="-terms" && i+1<argc) nterms=(long)std::stod(argv[i+1]);
		}
	std::cerr << "nthreads=" << nthreads << " nterms=" << nterms << std::endl;
	std::vector<std::thread> threads(nthreads);
	std::vector<datum> data(nthreads);
	for (int i=0;i<nthreads;i++) {
		data[i].start = i*nterms/nthreads;
		data[i].stop = (i+1)*nterms/nthreads;
		}
	for (int i=0;i<nthreads;i++) {
		threads[i]=std::thread(harm, std::ref(data[i]));
		}
	for(std::thread& thread : threads) thread.join();
	double sum=0;
	for (datum d : data) sum+=d.sum;
	std::cout << "sum=" << sum << std::endl;
}