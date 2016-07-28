#include <boost/dynamic_bitset.hpp>
#include <iostream>

int main(int, char*[]){
  boost::dynamic_bitset<> x(5);
  x[0] = 1;
  x[1] = 1;
  x[4] = 1;

  boost::dynamic_bitset<> y(6);
  y[0] = 1;
  y[1] = 1;
  y[2] = 1;

  for (boost::dynamic_bitset<>::size_type i = 0; i < x.size(); ++i)
  std::cout << x[i];
  std::cout << "\n";

  std::cout << "x=" << x << "\n";
  std::cout << "y=" << y << "\n";
  if(y.to_ulong() < x.to_ulong()){std::cout << "y < x\n";}else{std::cout << "!(y < x)\n";}

  std::string seq("ACTG");
  boost::dynamic_bitset<> z(2*seq.length());

  for(char &c : seq){
    z <<= 2;
    switch(c) {
      case 'A': break;
      case 'C': z[0] = 1;
      break;
      case 'T': z[1] = 1;
      break;
      case 'G': z[0] = 1; z[1] = 1;
      break;
    }
  }

  std::cout << "seq=" << seq << " lenght=" << seq.length() << " expected bin=" << "00011011" << " actual_bin=" << z << "\n";

  return EXIT_SUCCESS;


}
