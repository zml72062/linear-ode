#include "symsolver.h"


#define NORMAL(mat) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(mat)).normal())


symsolver_analytic::symsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, int order)
    : var(x), coeff(order + 1), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    for (int i = 0; i <= order; i++) {
        coeff[i] = GiNaC::matrix(N, N);
        sol[i] = GiNaC::matrix(N, N);
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bool exist = false;
            auto current_denom = A(i, j).normal().denom();
            for (auto& existing_denom: rho_eqs) {
                if (GiNaC::is_a<GiNaC::numeric>((existing_denom / current_denom).normal())) {
                    exist = true;
                    break;
                }
            }
            if (!exist)
                rho_eqs.append(current_denom);

            auto A_series = A(i, j).series(var, order + 1);
            for (int k = 0; k <= order; k++)
                coeff[k](i, j) = A_series.coeff(var, k).normal();
        }
    }

    for (int i = 0; i < N; i++)
        sol[0](i, i) = 1;
}


symsolver_analytic::symsolver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order)
    : var(x), coeff(order + 1), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    if ((int)Y0.rows() != N)
        throw std::runtime_error("Y0 and A have inconsistent shapes!");

    for (int i = 0; i <= order; i++) {
        coeff[i] = GiNaC::matrix(N, N);
        sol[i] = GiNaC::matrix(N, Y0.cols());
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bool exist = false;
            auto current_denom = A(i, j).normal().denom();
            for (auto& existing_denom: rho_eqs) {
                if (GiNaC::is_a<GiNaC::numeric>((existing_denom / current_denom).normal())) {
                    exist = true;
                    break;
                }
            }
            if (!exist)
                rho_eqs.append(current_denom);

            auto A_series = A(i, j).series(var, order + 1);
            for (int k = 0; k <= order; k++)
                coeff[k](i, j) = A_series.coeff(var, k).normal();
        }
    }

    int C = Y0.cols();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < C; j++)
            sol[0](i, j) = Y0(i, j);
}


GiNaC::matrix symsolver_analytic::solution() {
    int N = sol[0].rows(), C = sol[0].cols();
    int order = coeff.size() - 1;
    for (int k = 0; k < order; k++) {
        sol[k + 1] = GiNaC::matrix(N, C);
        for (int j = 0; j <= k; j++)
            sol[k + 1] = NORMAL(sol[k + 1].add(NORMAL(coeff[k - j].mul(sol[j]))));
        sol[k + 1] = NORMAL(sol[k + 1].mul_scalar(GiNaC::ex(1) / (k + 1)));
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


GiNaC::lst symsolver_analytic::radius_eqs() {
    return rho_eqs;
}


symsolver_regular::symsolver_regular(const GiNaC::matrix& A, const GiNaC::ex& x, int order, sym_regular_t& rstruct)
    : var(x), prstruct(&rstruct), coeff(order + 1), solutions(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    for (int i = 0; i <= order; i++) {
        coeff[i] = GiNaC::matrix(N, N);
        solutions[i] = GiNaC::matrix(N, N);
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bool exist = false;
            auto current_denom = A(i, j).normal().denom();
            for (auto& existing_denom: rho_eqs) {
                if (GiNaC::is_a<GiNaC::numeric>((existing_denom / current_denom).normal())) {
                    exist = true;
                    break;
                }
            }
            if (!exist)
                rho_eqs.append(current_denom);

            auto A_series = A(i, j).series(var, order + 1);
            for (int k = 0; k < order; k++)
                coeff[k + 1](i, j) = A_series.coeff(var, k).normal();

            coeff[0](i, j) = prstruct->J(i, j);
        }
    }
}


void symsolver_regular::fill_in_solution(int psol, int puppersol) {
    int order = coeff.size() - 1, N = solutions[0].rows();
    
    for (int i = 0; i < N; i++)
        solutions[0](i, psol) = (i == psol ? 1 : 0);
    
    GiNaC::ex diag = prstruct->J(psol, psol);

    for (int i = 1; i <= order; i++) {
        diag = NORMAL(diag + 1);

        GiNaC::matrix temp(N, 1), tempY(N, 1);
        if (puppersol != -1)
            for (int j = 0; j < N; j++)
                temp(j, 0) = -solutions[i](j, puppersol);
        
        for (int l = 1; l <= i; l++) {
            for (int j = 0; j < N; j++)
                tempY(j, 0) = solutions[i - l](j, psol);
            temp = NORMAL(temp.add(NORMAL(coeff[l].mul(tempY))));
        }

        GiNaC::matrix to_inv = NORMAL(GiNaC::ex_to<GiNaC::matrix>(GiNaC::unit_matrix(N)).mul_scalar(diag).sub(coeff[0])),
                      inv = NORMAL(to_inv.inverse());
        temp = NORMAL(inv.mul(temp));
        for (int j = 0; j < N; j++)
            solutions[i](j, psol) = temp(j, 0);
    }
}


GiNaC::matrix symsolver_regular::solution() {
    int order = coeff.size() - 1, N = solutions[0].rows();
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


GiNaC::lst symsolver_regular::radius_eqs() {
    return rho_eqs;
}

