#include "symdiffeq.h"


int main() {
    // Bessel equation of arbitrary order
    auto eq1 = symdiffeq::from_string("{{-1/x, (eps^2-x^2)/(x^2)},{1,0}}");
    auto res1 = eq1.solve(0, 10);
    std::cout << "solution matrix: " << res1.first << std::endl;
    std::cout << "convergence radius: " << res1.second << std::endl;
    
    // Legendre equation of arbitrary order
    auto eq2 = symdiffeq::from_string("{{2*x/(1-x^2), -eps*(eps+1)/(1-x^2)},{1,0}}");
    auto res2 = eq2.solve(0, 10);
    std::cout << "solution matrix: " << res2.first << std::endl;
    std::cout << "convergence radius: " << res2.second << std::endl;

    // resonant test case
    auto eq3 = symdiffeq::from_string("{{2*eps/(x^2-1),(1+eps)/(x-1)-x-1},{-1,(1+eps)/(x-1)}}");
    auto res3 = eq3.solve(1, 10);
    std::cout << "solution matrix: " << res3.first << std::endl;
    std::cout << "convergence radius: " << res3.second << std::endl;

    // special solution test case
    auto eq4 = symdiffeq::from_string("{{x/(2-(x-1-eps)^2), (2*x-1)/(x-1)},{x+1/(x+eps),-eps}}");
    auto res4 = eq4.solve(eq4.eps, GiNaC::matrix{{1-eq4.eps}, {0}}, 10);
    std::cout << "solution matrix: " << res4.first << std::endl;
    std::cout << "convergence radius: " << res4.second << std::endl;
}
