#ifndef SOLVER_H
#define SOLVER_H


#include "diffeq.h"


class solver_analytic {
public:
    /**
     * Get general solutions of [dY/dx = A(x) Y(x)] around origin, which is
     * assumed to be an analytic point of A(x).
     * 
     * @param A N by N matrix of coefficients
     * @param x differentiation variable
     * @param order series expansion order
     */
    solver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, int order);

    /**
     * Get a special solution of [dY/dx = A(x) Y(x)] around origin, which is
     * assumed to be an analytic point of A(x), with initial value [Y(0) = Y0].
     * 
     * @param A N by N matrix of coefficients
     * @param x differentiation variable
     * @param Y0 N by 1 matrix of initial value
     * @param order series expansion order
     */
    solver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order);

    /**
     * Return the solution around origin. If an initial value Y0 is given, 
     * will return an N by 1 matrix representing the special solution. 
     * Otherwise, will return an N by N matrix representing N basis vectors 
     * spanning the space of general solutions.
     * 
     * @param dig number of digits to keep
     */
    GiNaC::matrix solution(unsigned dig);

    /**
     * Return the radius of convergence.
     */
    GiNaC::ex radius();
private:
    GiNaC::ex     var;
    GiNaC::ex     rho;
    flint::ca_ctx ctx;
    std::vector<flint::ca_matrix> coeff;
    std::vector<flint::ca_matrix> sol;
};


class solver_regular {
public:
    /**
     * Get general solutions of [dY/dx = A(x) Y(x)] around origin, where we
     * assume that A(x) has an expansion of 
     * 
     *      A(x) = (Lambda / x) + A_analytic(x)
     *  
     * around origin, and "Lambda" is a Jordan canonical matrix with no 
     * eigenvalues that differ by a nonzero integer.
     * 
     * @param A N by N matrix, representing analytic part of the coefficient
     * A_analytic(x)
     * @param x differentiation variable
     * @param order series expansion order
     * @param rstruct reference to a `regular_t` struct, which saves singular
     * part of the coefficient
     */
    solver_regular(const GiNaC::matrix& A, const GiNaC::ex& x, int order, regular_t& rstruct);

    /**
     * Return the solution around origin. The solution is represented as an
     * N by N matrix, whose N columns represent N basis vectors spanning the
     * space of general solutions.
     * 
     * @param dig number of digits to keep
     */
    GiNaC::matrix solution(unsigned dig);

    /**
     * Return the radius of convergence.
     */
    GiNaC::ex radius();
private:
    GiNaC::ex     var;
    GiNaC::ex     rho;
    regular_t*    prstruct;
    std::vector<flint::ca_matrix> coeff;
    std::vector<flint::ca_matrix> solutions;

    void fill_in_solution(int psol, int puppersol);
};


#endif // SOLVER_H

