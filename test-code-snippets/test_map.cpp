// Some test code to experiment with the STL "map" class

#include <map>
#include <boost/dynamic_bitset.hpp>
#include <iostream>

int main(int, char*[]){

  boost::dynamic_bitset<> x(5);
  x[0] = 1;
  x[1] = 1;
  x[4] = 1;

  std::map <boost::dynamic_bitset<>, std::string> testmap;
  std::cout << "Size of map is:" << testmap.size() << "\n";

  // Insert one million items:
  for(unsigned long int i=200000000; i < 220000000; i++){
    boost::dynamic_bitset<> y(128, i);
    if(testmap.find(y) == testmap.end())
      testmap[y] = "String #" + std::to_string(i);
    else
      std::cout << "Found " << y << " more than once.\n";
  }

  std::cout << "Size of map is:" << testmap.size() << "\n";

}
