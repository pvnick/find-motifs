#include <iostream>
#include <cmath>

int main() {
unsigned long long target = pow(1600000, 2);
std::cout << target << std::endl;
for (unsigned long long i = 0; i < target; ++i) {
if ((i%10000000) == 0) {
std::cout << i << std::endl;
}
}
return 0;
}
