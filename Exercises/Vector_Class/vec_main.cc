#include <iostream>
#include "vec_class.h"

int main(){
    vector3D vec1 = vector3D(3, 4, 0);
    vector3D vec2 = vector3D(0, 7, 12);
    std::cout << (2*(vec1 - vec2)).z << std::endl;
    vec2.normalize();
    std::cout << vec2.z << " should equal " << 12/std::sqrt(12*12 + 7*7) << std::endl;
    return 0;
};