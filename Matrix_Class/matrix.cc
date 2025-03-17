#include"matrix.h"
#include<string>
#include<algorithm>
#include<cmath>
#define SELF (*this)
#define FORV(i,v) for(int i=0;i<v.size();i++)
#define FOR_COLS(i,A) for(int i=0;i<A.size2();i++)
namespace pp{

bool approx(NUMBER x,NUMBER y,NUMBER acc=1e-6,NUMBER eps=1e-6){
	if(std::fabs(x-y) < acc)return true;
	if(std::fabs(x-y) < eps*(std::fabs(x)+std::fabs(y)))return true;
	return false;
}

bool approx(const vector& u,const vector& v,NUMBER acc,NUMBER eps){
	if(u.size()!=v.size())return false;
	for(int i=0;i<u.size();i++)if(!approx(u[i],v[i],acc,eps))return false;
	return true;
}

vector& vector::operator+=(const vector& other) {
	FORV(i,SELF) data[i]+=other.data[i];
	return SELF; }

vector& vector::operator-=(const vector& other) {
	FORV(i,SELF) data[i]-=other.data[i];
	return SELF; }

vector& vector::operator*=(NUMBER x) {
	FORV(i,SELF) data[i]*=x;
	return SELF; }

vector& vector::operator/=(NUMBER x) {
	FORV(i,SELF) data[i]/=x;
	return SELF; }

void vector::print(std::string s,FILE* stream) const {
	fprintf(stream,"%s\n",s.c_str());
	FORV(i,SELF)fprintf(stream,"%9.4g ",(double)SELF[i]);
	fprintf(stream,"\n");
	}

vector operator/(const vector& v, NUMBER x){
	vector r=v;
	r/=x;
	return r; }

vector operator*(const vector& v, NUMBER x){
	vector r=v;
	r*=x;
	return r; }

vector operator*(NUMBER x,const vector& a){ return a*x; }

vector operator+(const vector& a, const vector& b){
	vector r=a;
	r+=b;
	return r; }

vector operator-(const vector& a, const vector& b){
	vector r=a;
	r-=b;
	return r; }

void matrix::resize(int n, int m){
	cols.resize(m);
	for(int i=0;i<m;++i)cols[i].resize(n);
	}

matrix matrix::transpose(){
    //matrix R; R.resize(size2(),size1());
	matrix R(size2(),size1());
    for(int j=0;j<R.size2();++j)
    for(int i=0;i<R.size1();++i)
        R(i,j)=SELF(j, i);
    return R;
    }

matrix& matrix::operator+=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]+=other[i];
	return SELF; }

matrix& matrix::operator-=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]-=other[i];
	return SELF; }

matrix& matrix::operator*=(NUMBER x) {
	FOR_COLS(i,SELF) SELF[i]*=x;
	return SELF; }

matrix& matrix::operator/=(NUMBER x) {
	FOR_COLS(i,SELF) SELF[i]/=x;
	return SELF; }

matrix operator/(const matrix& A,NUMBER x){
	matrix R=A;
	R/=x;
	return R; }

matrix operator*(const matrix& A,NUMBER x){
	matrix R=A;
	R*=x;
	return R; }

matrix operator*(NUMBER x,const matrix& A){
	return A*x; }

matrix operator+(const matrix& A, const matrix& B){
	matrix R=A;
	R+=B;
	return R; }

matrix operator-(const matrix& A, const matrix& B){
	matrix R=A;
	R-=B;
	return R; }

vector operator*(const matrix& M, const vector& v){
	vector r; r.resize(M.size1());
	for(int i=0;i<r.size();i++){
		NUMBER sum=0;
		for(int j=0;j<v.size();j++)sum+=M(i, j)*v[j];
		r[i]=sum;
		}
	return r;
	}

matrix operator*(const matrix& A, const matrix& B){
	matrix R(A.size1(),B.size2());
	for(int k=0;k<A.size2();k++)
	for(int j=0;j<B.size2();j++)
		{
		for(int i=0;i<A.size1();i++)R(i, j)+=A(i,k)*B(k,j);
		}
	return R;
	}
/*
vector matrix::get_col(int j){
	vector cj=SELF[j];
	return cj;
	}

void matrix::set_col(int j,vector& cj){
	SELF[i]=cj;
	}
*/

void matrix::print(std::string s,FILE* stream){
	fprintf(stream,"%s\n",s.c_str());
	for(int i=0;i<size1();i++){
		for(int j=0;j<size2();j++)fprintf(stream,"%9.4g ",(double)SELF(i, j));
		fprintf(stream,"\n");
		}
	}
}//pp