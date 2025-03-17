#include"matrix.h"
pp::matrix eye(int n){
	pp::matrix M(n,n);
	for(int i=0;i<n;i++)M(i,i)=1;
	return M;
	}
int main(){
	int n=3;
	pp::vector v(n);
	v.print("v=");
	for(int i=0;i<v.size();i++)v[i]=i+1;
	v.print("v=");
	pp::vector u=v;
	u.print("u=");
	v[0]=999;
	u.print("u=");
	pp::vector w(v);
	w.print("w=");
	w+=v;
	w.print("w+=v");
	(w/10).print("w/10");
	(u+u).print("u+u");
	pp::matrix M(n,n);
	M.print("M=");
	for(int i=0;i<M.size1();i++)
	for(int j=0;j<M.size2();j++)
		M(i,j)=i+j;
	M.print("M=");
	(M/10).print("M/10");
	(M+M).print("M+M");
	(M*(M+M)*M).print("M*(M+M)*M");
	pp::matrix O=eye(n);
	O.print("O=");
return 0;
}