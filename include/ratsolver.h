#ifndef RATSOLVER_H
#define RATSOLVER_H


/**
 * ratsolver.h - Optimized solver for rational coefficients.
 */


#include "symdiffeq.h"


class ratsolver_analytic {
public:
    /**
     * Get general solutions of [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] around 
     * origin, which is assumed to be an analytic point of A(x, eps).
     * 
     * @param A N by N matrix of coefficients
     * @param x differentiation variable
     * @param order series expansion order
     */
    ratsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, int order);

    /**
     * Get a special solution of [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] around 
     * origin, which is assumed to be an analytic point of A(x, eps), with 
     * initial value [Y(0, eps) = Y0(eps)].
     * 
     * @param A N by N matrix of coefficients
     * @param x differentiation variable
     * @param Y0 N by 1 matrix of initial value
     * @param order series expansion order
     */
    ratsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order);

    /**
     * Return the solution around origin. If an initial value Y0(eps) is given, 
     * will return an N by 1 matrix representing the special solution. 
     * Otherwise, will return an N by N matrix representing N basis vectors 
     * spanning the space of general solutions.
     */
    GiNaC::matrix solution();

    /**
     * Return an algebraic equation whose roots mark the singularities of 
     * the equation.
     */
    GiNaC::ex radius_eq();
private:
    GiNaC::ex     var;
    GiNaC::ex     rho_eq;
    std::vector<GiNaC::ex>     denom;
    std::vector<GiNaC::matrix> numer;
    std::vector<GiNaC::matrix> sol;
};


class ratsolver_regular {
public:
    /**
     * Get general solutions of [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] around 
     * origin, where we assume that A(x, eps) has an expansion of 
     * 
     *      A(x, eps) = (Lambda(eps) / x) + A_analytic(x, eps)
     *  
     * around origin, and "Lambda(eps)" is a Jordan canonical matrix with no 
     * eigenvalues that differ by a nonzero integer.
     * 
     * @param A N by N matrix, representing analytic part of the coefficient
     * A_analytic(x, eps)
     * @param x differentiation variable
     * @param order series expansion order
     * @param rstruct reference to a `sym_regular_t` struct, which saves 
     * singular part of the coefficient
     */
    ratsolver_regular(const GiNaC::matrix& A, const GiNaC::ex& x, int order, sym_regular_t& rstruct);

    /**
     * Return the solution around origin. The solution is represented as an
     * N by N matrix, whose N columns represent N basis vectors spanning the
     * space of general solutions.
     */
    GiNaC::matrix solution();

    /**
     * Return an algebraic equation whose roots mark the singularities of 
     * the equation.
     */
    GiNaC::ex radius_eq();
private:
    GiNaC::ex     var;
    GiNaC::ex     rho_eq;
    sym_regular_t* prstruct;
    std::vector<GiNaC::ex>     denom;
    std::vector<GiNaC::matrix> numer;
    std::vector<GiNaC::matrix> solutions;

    void fill_in_solution(int psol, int puppersol);
};


#endif // RATSOLVER_H

