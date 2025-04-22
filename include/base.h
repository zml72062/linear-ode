#ifndef BASE_H
#define BASE_H


#include <ginac/ginac.h>


typedef struct moser_struct_t {
    int negative_power;
    GiNaC::matrix A0;
    GiNaC::matrix A1;

    moser_struct_t(int rows): A0(rows, rows), A1(rows, rows) { }

    GiNaC::ex moser_rank() {
        if (negative_power <= 1)
            return 1;
        return negative_power - 1 + GiNaC::ex(A0.rank()) / A0.rows();
    }
} moser_t;


class base_diffeq {
public:
    base_diffeq(): x("x") { }
    void initialize(const GiNaC::matrix& _matrix);
    
    int N() { return raw_coeff.rows(); }
    const GiNaC::matrix& get_raw_coeff() { return raw_coeff; }
    const GiNaC::matrix& get_coeff()     { return coeff;     }
    const GiNaC::matrix& get_transform() { return transform; }

    moser_t get_moser_struct(const GiNaC::ex& x0);
    void moser_reduction_one_step(const GiNaC::ex& x0);
    void moser_reduction(const GiNaC::ex& x0);

    bool is_analytic(const GiNaC::ex& x0);
    bool is_first_order(const GiNaC::ex& x0);
    bool is_regular(const GiNaC::ex& x0);

    void update(const GiNaC::matrix& new_transform);
    void clear_all_reduction();
protected:
    GiNaC::symbol x;
    GiNaC::matrix raw_coeff;
    GiNaC::matrix coeff;
    GiNaC::matrix transform;
};


#endif // BASE_H

