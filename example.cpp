#include "diffeq.h"


int main() {
    // Bessel equation of order 3/2
    GiNaC::Digits = 100;
    auto eq1 = diffeq::from_string("{{-1/x, (9/4-x^2)/(x^2)},{1,0}}");
    auto res1 = eq1.solve(0, 10, 100);
    std::cout << "solution matrix: " << res1.first << std::endl;
    std::cout << "convergence radius: " << res1.second << std::endl;

    // Bessel equation of order 2, centered at x=-1
    auto eq2 = diffeq::from_string("{{-1/(x+1), (4-(x+1)^2)/((x+1)^2)},{1,0}}");
    auto res2 = eq2.solve(-1, 10, 100);
    std::cout << "solution matrix: " << res2.first << std::endl;
    std::cout << "convergence radius: " << res2.second << std::endl;

    // Legendre equation of order 2
    auto eq3 = diffeq::from_string("{{2*x/(1-x^2), -6/(1-x^2)},{1,0}}");
    auto res3 = eq3.solve(0, 10, 100);
    std::cout << "solution matrix: " << res3.first << std::endl;
    std::cout << "convergence radius: " << res3.second << std::endl;
}
