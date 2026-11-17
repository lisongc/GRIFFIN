// Compile ONLY against independently extracted upstream v1.1.1 sources.
// Plain inval receives already resolved masses/widths: no second EW solve.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "EWPOZ2.h"
#include "xscaas.h"
using namespace griffin;

int main()
{
  try {
    const int indices[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,20,21,22,23};
    std::cout << std::setprecision(17);
    std::string line;
    while(std::getline(std::cin,line))
    {
      if(line.empty()) continue;
      std::stringstream stream(line);
      std::string value;
      std::vector<std::string> row;
      while(std::getline(stream,value,',')) row.push_back(value);
      if(row.size()!=26) throw std::invalid_argument("expected 26 reference fields");
      if(row[0].empty() || row[0].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
        throw std::invalid_argument("invalid id");
      const int final=std::stoi(row[1]);
      const double s=std::stod(row[2]), c=std::stod(row[3]);
      inval input;
      for(unsigned k=0;k<sizeof(indices)/sizeof(int);++k)
        input.set(indices[k],std::stod(row[4+k]));
      FA_SMNNLO fai(ELE,input), faf(final,input);
      SW_SMNNLO swi(ELE,input), swf(final,input);
      mat_SMaas matrix(ELE,final,VEC,VEC,fai,faf,swi,swf,s,c,input);
      const int left[]={VEC,AXV,VEC,AXV}, right[]={VEC,VEC,AXV,AXV};
      std::array<Cplx,4> a;
      for(unsigned k=0;k<4;++k)
      {
        matrix.setform(left[k],right[k]); a[k]=matrix.result();
        if(!std::isfinite(real(a[k])) || !std::isfinite(imag(a[k])))
          throw std::runtime_error("nonfinite upstream amplitude: "+row[0]);
      }
      const double nc=final<=5 ? 3.0 : 1.0;
      // Independent copy of the original upstream examples' angular square.
      const double angular=real((1+c*c)*(a[0]*conj(a[0])+a[1]*conj(a[1])
          +a[2]*conj(a[2])+a[3]*conj(a[3]))
          +4*c*(a[0]*conj(a[3])+a[2]*conj(a[1])));
      const double sigma=nc*s/(32*Pi)*0.38937966e9*angular;
      if(!std::isfinite(sigma)) throw std::runtime_error("nonfinite cross section");
      std::cout << "{\"id\":\"" << row[0] << "\",\"amplitudes\":[";
      for(unsigned k=0;k<4;++k)
      {
        if(k) std::cout << ',';
        std::cout << '[' << real(a[k]) << ',' << imag(a[k]) << ']';
      }
      std::cout << "],\"sigma_pb\":" << sigma << "}\n";
    }
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "independent upstream reference: " << error.what() << '\n';
    return 1;
  }
}
