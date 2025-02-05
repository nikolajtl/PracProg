#include <iostream>
#include <vector>

int main(){
    double x=5;
    double y=x;
    double& z=x;
    int n=3;
    std::vector<double> a(3);
    for (int i=0; i<a.size(); i++){
        a[i]=i+1;
    }
    std::vector<double>& b=a;
    std::cout << b << std::endl;
    return 0;
}
