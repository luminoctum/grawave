#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <vector>

using namespace std;
using namespace Eigen;

int main(){
    int i = 1;
    int j = 3;
    //ArrayXXf a(1,3); // ERROR
    //MatrixXf a(1,3); // FFT fail
    Array<float,1,3> a;
    //Matrix<float,1,3> a;
    //Matrix<complex<float>,1,3> b;
    Array<complex<float>,1,3> b;
    //ArrayXXcf b(1,3);
    //MatrixXcf b(1,3);
    a << 1,2,3;
    b << 1,2,3;
    FFT<float> fft;
    fft.fwd(b, a.matrix());
    cout << a << endl;
    cout << b.array() << endl;
}
