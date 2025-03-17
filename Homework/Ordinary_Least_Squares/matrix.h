#ifndef HAVE_MATRIX_H
#define HAVE_MATRIX_H
#ifdef LONG_DOUBLE
	#define NUMBER long double
#else
	#define NUMBER double
#endif
#include<string>
#include<vector>
namespace pp{
struct vector {
	std::vector<NUMBER> data;
	vector(int n) : data(n) {}
	vector()			=default;
	vector(const vector&)		=default;
	vector(vector&&)		=default;
	vector(const std::vector<NUMBER>& v) : data(v) {}
	vector& operator=(const vector&)=default;
	vector& operator=(vector&&)	=default;
	int size() const {return data.size();}
	void resize(int n) {data.resize(n);}
	void push_back(const NUMBER& val) {data.push_back(val); } // Add an entry to the vector
	NUMBER norm() const;
	NUMBER dot(const vector& other) const;
	NUMBER& operator[](int i) {return data[i];}
	const NUMBER& operator[](int i) const {return data[i];}
	vector& operator+=(const vector&);
	vector& operator-=(const vector&);
	vector& operator*=(NUMBER);
	vector& operator/=(NUMBER);
	bool operator==(const vector& other) const; // Vector equality operator
	void print(std::string s="",FILE* stream=stdout) const;
};

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const vector&, NUMBER);
vector operator*(NUMBER, const vector&);
vector operator/(const vector&, NUMBER);
bool approx(const vector&, const vector&, NUMBER acc=1e-6, NUMBER eps=1e-6);

struct matrix {
	std::vector<vector> cols;
	matrix()=default;
	matrix(int n,int m) : cols(m, vector(n)) {}
	matrix(const matrix& other)=default;
	matrix(matrix&& other)=default;
	matrix(const std::vector<vector>& v) : cols(v) {} // Constructor via. feeding column vectors directly
	matrix& operator=(const matrix& other)=default;
	matrix& operator=(matrix&& other)=default;
	int size1() const {return cols.empty() ? 0 : cols[0].size(); }
	int size2() const {return cols.size();}
	void resize(int n, int m);
	void push_back(const vector& col) { cols.push_back(col); } // Add column vectors by pushing back
	
	NUMBER get (int i, int j) {return cols[j][i];}
	void set(int i, int j, NUMBER value){cols[j][i] = value;}
	NUMBER& operator()(int i, int j){return cols[j][i];}
	const NUMBER& operator()(int i, int j) const {return cols[j][i];}
	// NUMBER& operator[](int i, int j){return cols[j][i];}
	// const NUMBER& operator[](int i, int j) const {return cols[j][i];}
	vector& operator[](int i){return cols[i];}
	const vector& operator[](int i) const {return cols[i];}
//	vector get_col(int j);
//	void set_col(int j,vector& cj);
	matrix transpose();

	matrix& operator+=(const matrix&);
	matrix& operator-=(const matrix&);
	matrix& operator*=(const matrix&);
	matrix& operator*=(const NUMBER);
	matrix& operator/=(const NUMBER);
	matrix  operator^(int);
	bool operator==(const matrix& other) const; // Matrix equality operator

	//void print(const char* s="");
	void print(std::string s="",FILE* stream=stdout);
};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, NUMBER);
matrix operator*(NUMBER, const matrix&);
matrix operator/(const matrix&, NUMBER);
vector operator*(const matrix&, const vector&);
}
#endif