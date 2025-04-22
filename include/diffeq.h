#ifndef DIFFEQ_H
#define DIFFEQ_H


#include "wrapper.h"
#include "base.h"


typedef struct regular_struct_t {
    flint::ca_ctx    ctx;
    flint::ca_matrix A0;
    flint::ca_matrix J;
    flint::ca_matrix P;
    flint::ca_vec    lambdas;
    std::vector<ulong> mult;
    slong            num_blocks;
    std::vector<slong> block_lambda;
    std::vector<slong> block_size;
    std::vector<std::vector<std::pair<int, int>>> eiggroups;
    std::vector<int> block_structure;
    std::vector<int> block_structure_idx;

    void initialize(const GiNaC::matrix& _A0);
    void group_eigenvalues();
    void update();
} regular_t;


class diffeq: public base_diffeq {
public:
    diffeq() { }
    
    /**
     * Define a homogeneous linear differential equation 
     * [dY/dx = A(x) Y(x)] by a square matrix A(x).
     * 
     * @param _coeff a square matrix A(x) of coefficients
     * @param _var differentiation variable x
     */
    diffeq(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var);

    /**
     * Define a homogeneous linear differential equation 
     * [dY/dx = A(x) Y(x)] from the string representation of A(x).
     * 
     * @param _str string representation of coefficient matrix A(x)
     * in Mathematica form
     */
    static diffeq from_string(const std::string& _str);
    
    int regular_reduction_one_step(const GiNaC::ex& x0, unsigned digits);
    void regular_reduction(const GiNaC::ex& x0, unsigned digits);

    /**
     * Give the general solution to the differential equation 
     * [dY/dx = A(x) Y(x)] around x=x0 by series expansion method, 
     * assuming that x0 is an analytic point or a regular singular point. 
     * 
     * @param x0 point of evaluation, assumed to be an analytic point or
     * a regular singular point of the equation
     * @param order series expansion order
     * @param digits number of digits to keep for the coefficients
     * 
     * @returns a pair consisting of
     *      (a) a square matrix having the same shape as A(x), each column
     *          representing a linearly independent solution
     *      (b) the radius of convergence
     * 
     * @throws throws "not implemented error", if x0 is neither an analytic 
     * point nor a regular singular point (a.k.a. an irregular singular point)
     */
    std::pair<GiNaC::matrix, GiNaC::ex> solve(const GiNaC::ex& x0, unsigned order, unsigned digits);

    /**
     * Give the special solution to the differential equation
     * [dY/dx = A(x) Y(x)] with initial condition [Y(x0) = Y0] around x=x0
     * by series expansion method, assuming that x0 is an analytic point.
     * 
     * @param x0 point of evaluation, assumed to be an analytic point of the
     * equation
     * @param Y0 the initial value Y(x0)
     * @param order series expansion order
     * @param digits number of digits to keep for the coefficients
     * 
     * @returns a pair consisting of 
     *      (a) a matrix having the same shape as Y0, representing the special
     *          solution to the initial value problem
     *      (b) the radius of convergence
     * 
     * @throws throws "not implemented error", if x0 is not an analytic point
     * (a.k.a. a singularity)
     */
    std::pair<GiNaC::matrix, GiNaC::ex> solve(const GiNaC::ex& x0, const GiNaC::matrix& Y0, unsigned order, unsigned digits);
private:
    // only for regular singularity
    regular_t     reg_struct;
};



#endif // DIFFEQ_H

