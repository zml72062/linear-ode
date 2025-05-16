#include "ratsolver.h"


#define NORMAL(mat) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(mat)).normal())


ratsolver_analytic::ratsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, int order)
    : var(x), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    
    GiNaC::ex common_denom = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto current_denom = A(i, j).normal().denom();
            common_denom = GiNaC::lcm(common_denom, current_denom);
        }
    }

    GiNaC::ex expanded_common_denom = common_denom.expand().collect(var);
    if (expanded_common_denom.ldegree(var) > 0)
        throw std::runtime_error("x=0 is not an analytic point of A");
    
    rho_eq = common_denom;
    int d_degree = expanded_common_denom.degree(var);
    denom = std::vector<GiNaC::ex>(d_degree + 1);
    for (int i = 0; i <= d_degree; i++)
        denom[i] = expanded_common_denom.coeff(var, i);

    GiNaC::matrix numerator_polys(N, N);
    int n_degree = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            numerator_polys(i, j) = (A(i, j) * common_denom).normal().expand().collect(var);
            int current_n_degree = numerator_polys(i, j).degree(var);
            if (current_n_degree > n_degree)
                n_degree = current_n_degree;
        }
    }

    numer = std::vector<GiNaC::matrix>(n_degree + 1, GiNaC::matrix(N, N));
    for (int d = 0; d <= n_degree; d++)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                numer[d](i, j) = numerator_polys(i, j).coeff(var, d);

    for (int i = 0; i <= order; i++)
        sol[i] = GiNaC::matrix(N, N);
    
    for (int i = 0; i < N; i++)
        sol[0](i, i) = 1;
}


ratsolver_analytic::ratsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order)
    : var(x), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    if ((int)Y0.rows() != N)
        throw std::runtime_error("Y0 and A have inconsistent shapes!");
    
    GiNaC::ex common_denom = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto current_denom = A(i, j).normal().denom();
            common_denom = GiNaC::lcm(common_denom, current_denom);
        }
    }

    GiNaC::ex expanded_common_denom = common_denom.expand().collect(var);
    if (expanded_common_denom.ldegree(var) > 0)
        throw std::runtime_error("x=0 is not an analytic point of A");

    rho_eq = common_denom;
    int d_degree = expanded_common_denom.degree(var);
    denom = std::vector<GiNaC::ex>(d_degree + 1);
    for (int i = 0; i <= d_degree; i++)
        denom[i] = expanded_common_denom.coeff(var, i);

    GiNaC::matrix numerator_polys(N, N);
    int n_degree = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            numerator_polys(i, j) = (A(i, j) * common_denom).normal().expand().collect(var);
            int current_n_degree = numerator_polys(i, j).degree(var);
            if (current_n_degree > n_degree)
                n_degree = current_n_degree;
        }
    }

    numer = std::vector<GiNaC::matrix>(n_degree + 1, GiNaC::matrix(N, N));
    for (int d = 0; d <= n_degree; d++)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                numer[d](i, j) = numerator_polys(i, j).coeff(var, d);

    for (int i = 0; i <= order; i++)
        sol[i] = GiNaC::matrix(N, Y0.cols());
    
    int C = Y0.cols();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < C; j++)
            sol[0](i, j) = Y0(i, j);
}


GiNaC::matrix ratsolver_analytic::solution() {
    int N = sol[0].rows(), C = sol[0].cols();
    int order = sol.size() - 1;
    int L_D = denom.size() - 1, L_N = numer.size() - 1;

    for (int k = 0; k < order; k++) {
        auto n_sum = GiNaC::matrix(N, C);
        int n_terms = std::min(k, L_N);
        for (int j = 0; j <= n_terms; j++)
            n_sum = NORMAL(n_sum.add(NORMAL(numer[j].mul(sol[k - j]))));
        
        auto d_sum = GiNaC::matrix(N, C);
        int d_terms = std::min(k, L_D);
        for (int j = 1; j <= d_terms; j++)
            d_sum = NORMAL(d_sum.add(NORMAL(sol[k - j + 1].mul_scalar(denom[j] * (k - j + 1)))));
        
        sol[k + 1] = NORMAL(n_sum.sub(d_sum).mul_scalar(GiNaC::ex(1) / denom[0] / (k + 1)));
    }

    GiNaC::matrix M(N, C);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < C; j++) {
            M(i, j) = 0;
            for (int k = 0; k <= order; k++) {
                M(i, j) += (sol[k](i, j) * GiNaC::pow(var, k));
            }
        }
    }

    return M;
}


GiNaC::ex ratsolver_analytic::radius_eq() {
    return rho_eq;
}


ratsolver_regular::ratsolver_regular(const GiNaC::matrix& A, const GiNaC::ex& x, int order, sym_regular_t& rstruct)
    : var(x), prstruct(&rstruct), solutions(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");

    GiNaC::ex common_denom = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto current_denom = A(i, j).normal().denom();
            common_denom = GiNaC::lcm(common_denom, current_denom);
        }
    }

    GiNaC::ex expanded_common_denom = common_denom.expand().collect(var);
    if (expanded_common_denom.ldegree(var) > 0)
        throw std::runtime_error("x=0 is not an analytic point of A");
    
    rho_eq = common_denom;
    int d_degree = expanded_common_denom.degree(var);
    denom = std::vector<GiNaC::ex>(d_degree + 1);
    for (int i = 0; i <= d_degree; i++)
        denom[i] = expanded_common_denom.coeff(var, i);

    GiNaC::matrix numerator_polys(N, N);
    int n_degree = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            numerator_polys(i, j) = (A(i, j) * common_denom).normal().expand().collect(var);
            int current_n_degree = numerator_polys(i, j).degree(var);
            if (current_n_degree > n_degree)
                n_degree = current_n_degree;
        }
    }

    numer = std::vector<GiNaC::matrix>(n_degree + 1, GiNaC::matrix(N, N));
    for (int d = 0; d <= n_degree; d++)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                numer[d](i, j) = numerator_polys(i, j).coeff(var, d);

    for (int i = 0; i <= order; i++)
        solutions[i] = GiNaC::matrix(N, N);
}


void ratsolver_regular::fill_in_solution(int psol, int puppersol) {
    int order = solutions.size() - 1, N = solutions[0].rows();
    int L_D = denom.size() - 1, L_N = numer.size() - 1;

    for (int i = 0; i < N; i++)
        solutions[0](i, psol) = (i == psol ? 1 : 0);
    
    GiNaC::ex diag = prstruct->J(psol, psol);

    for (int i = 1; i <= order; i++) {
        diag = NORMAL(diag + 1);

        GiNaC::matrix part0(N, 1), part1(N, 1), part2(N, 1);
        GiNaC::matrix temp1(N, 1), temp21(N, 1), temp22(N, 1);
        if (puppersol != -1)
            for (int j = 0; j < N; j++)
                part0(j, 0) = -solutions[i](j, puppersol);
        
        int terms1 = std::min(i - 1, L_N);
        for (int l = 0; l <= terms1; l++) {
            for (int j = 0; j < N; j++)
                temp1(j, 0) = solutions[i - l - 1](j, psol);
            part1 = NORMAL(part1.add(NORMAL(numer[l].mul(temp1))));
        }
        part1 = NORMAL(part1.mul_scalar(GiNaC::ex(1) / denom[0]));

        int terms2 = std::min(i, L_D);
        for (int l = 1; l <= terms2; l++) {
            if (puppersol != -1)
                for (int j = 0; j < N; j++)
                    temp21(j, 0) = solutions[i - l](j, puppersol);
            
            for (int j = 0; j < N; j++)
                temp22(j, 0) = solutions[i - l](j, psol);
            
            part2 = NORMAL(part2.add(
                NORMAL(temp21.add(
                    NORMAL(NORMAL(GiNaC::ex_to<GiNaC::matrix>(GiNaC::unit_matrix(N))
                        .mul_scalar(diag - l).sub(prstruct->J)
                    ).mul(temp22)))).mul_scalar(denom[l])
                )
            );
        }
        part2 = NORMAL(part2.mul_scalar(GiNaC::ex(1) / denom[0]));

        auto result = NORMAL(
            NORMAL(NORMAL(
                GiNaC::ex_to<GiNaC::matrix>(GiNaC::unit_matrix(N))
                .mul_scalar(diag).sub(prstruct->J)
            ).inverse()).mul(
                NORMAL(part0.add(part1.sub(part2)))
            )
        );

        for (int j = 0; j < N; j++)
            solutions[i](j, psol) = result(j, 0);
    }
}


GiNaC::matrix ratsolver_regular::solution() {
    int order = solutions.size() - 1, N = solutions[0].rows();
    fill_in_solution(0, -1);
    for (int i = 1; i < N; i++) {
        if(prstruct->J(i - 1, i) == 1)
            fill_in_solution(i, i - 1);
        else
            fill_in_solution(i, -1);
    }

    GiNaC::matrix sol(N, N), finalsol(N, N);
    int blockstart = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            sol(i, j) = 0;
            for (int k = 0; k <= order; k++)
                sol(i, j) += (solutions[k](i, j) * GiNaC::pow(var, k));
        }

        if (j == 0 || !(bool)(prstruct->J(j - 1, j) == 1)) {
            for (int i = 0; i < N; i++)
                finalsol(i, j) = sol(i, j) * GiNaC::pow(var, prstruct->J(j, j));
            blockstart = j;
        } else {
            for (int i = 0; i < N; i++) {
                finalsol(i, j) = sol(i, j);
                int k = j;
                while (k != blockstart) {
                    k--;
                    finalsol(i, j) += sol(i, k) * GiNaC::pow(GiNaC::log(var), j - k) / GiNaC::tgamma(j - k + 1);
                }
                finalsol(i, j) *= GiNaC::pow(var, prstruct->J(j, j));
            }
        }
    }
    return finalsol;
}


GiNaC::ex ratsolver_regular::radius_eq() {
    return rho_eq;
}

