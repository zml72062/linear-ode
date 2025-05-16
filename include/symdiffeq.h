#ifndef SYMDIFFEQ_H
#define SYMDIFFEQ_H


#include "base.h"


typedef struct symbolic_regular_struct_t {
    GiNaC::matrix A0;
    GiNaC::matrix J;
    GiNaC::matrix P;
    GiNaC::lst    lambdas;
    std::vector<int> mult;
    int           num_blocks;
    std::vector<int> block_lambda;
    std::vector<int> block_size;
    std::vector<std::vector<std::pair<int, int>>> eiggroups;
    std::vector<int> block_structure;
    std::vector<int> block_structure_idx;

    void initialize(const GiNaC::matrix& _A0);
    void group_eigenvalues();
    void update();
} sym_regular_t;


class symdiffeq: public base_diffeq {
public:
    symdiffeq(): eps("eps") { }

    /**
     * Define a symbolic homogeneous linear differential equation 
     * [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] by a square matrix 
     * A(x, eps).
     * 
     * Notice that solving a symbolic linear differential equation
     * can be significantly slower than solving a non-symbolic one
     * (i.e., one with no dependence on the parameter `eps`). Thus
     * it is advised that one instantiates a `diffeq` instead of a
     * `symdiffeq` when the coefficient matrix A(x) is independent
     * from `eps`.
     * 
     * @param _coeff a square matrix A(x, eps) of coefficients
     * @param _var differentiation variable x
     * @param _eps an additional parameter eps
     */
    symdiffeq(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps);

    /**
     * Define a symbolic homogeneous linear differential equation 
     * [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] from the string
     * representation of A(x, eps).
     * 
     * @param _str string representation of coefficient matrix 
     * A(x, eps) in Mathematica form
     */
    static symdiffeq from_string(const std::string& _str);

    int regular_reduction_one_step(const GiNaC::ex& x0);
    void regular_reduction(const GiNaC::ex& x0);

    /**
     * Give the general solution to the differential equation 
     * [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] around x=x0 by series 
     * expansion method, assuming that x0 is an analytic point or a 
     * regular singular point. 
     * 
     * @param x0 point of evaluation, assumed to be an analytic point or
     * a regular singular point of the equation
     * @param order series expansion order
     * @param opt whether to use an algorithm optimized for rational coefficient
     * matrices (default `true`)
     * 
     * @returns a pair consisting of
     *      (a) a square matrix having the same shape as A(x, eps), each 
     *          column representing a linearly independent solution
     *      (b) a list of polynomial equations with x being the variable,
     *          whose roots mark the singularities of the equation minus x0
     * 
     * @throws throws "not implemented error", if x0 is neither an analytic 
     * point nor a regular singular point (a.k.a. an irregular singular point)
     */
    std::pair<GiNaC::matrix, GiNaC::lst> solve(const GiNaC::ex& x0, unsigned order, bool opt = true);

    /**
     * Give the special solution to the differential equation
     * [d/dx Y(x, eps) = A(x, eps) Y(x, eps)] with initial condition 
     * [Y(x0, eps) = Y0(eps)] around x=x0 by series expansion method, 
     * assuming that x0 is an analytic point.
     * 
     * @param x0 point of evaluation, assumed to be an analytic point of the
     * equation
     * @param Y0 the initial value Y(x0, eps)
     * @param order series expansion order
     * @param opt whether to use an algorithm optimized for rational coefficient
     * matrices (default `true`)
     * 
     * @returns a pair consisting of 
     *      (a) a matrix having the same shape as Y0, representing the special
     *          solution to the initial value problem
     *      (b) a list of polynomial equations with x being the variable,
     *          whose roots mark the singularities of the equation minus x0
     * 
     * @throws throws "not implemented error", if x0 is not an analytic point
     * (a.k.a. a singularity)
     */
    std::pair<GiNaC::matrix, GiNaC::lst> solve(const GiNaC::ex& x0, const GiNaC::matrix& Y0, unsigned order, bool opt = true);
    
    GiNaC::symbol eps;
private:
    // only for regular singularity
    sym_regular_t reg_struct;
};



#endif // SYMDIFFEQ_H

