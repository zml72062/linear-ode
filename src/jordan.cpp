#include "jordan.h"
#include <vector>


#define SUBMATRIX(m, r, nr, c, nc) GiNaC::ex_to<GiNaC::matrix>(GiNaC::sub_matrix((m), (r), (nr), (c), (nc)))
#define NORMAL(mat) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(mat)).normal())
#define SUBS(m, r) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).subs((r), GiNaC::subs_options::algebraic))


/**
 * Add all factors of `poly` to `lst`.
 */
static void set_factors(const GiNaC::ex& poly, GiNaC::lst& lst) {
    if (!GiNaC::is_exactly_a<GiNaC::mul>(poly)) {
        lst.append(poly);
        return;
    }

    for (auto iter = poly.begin(); iter != poly.end(); ++iter)
        set_factors(*iter, lst);
}

static GiNaC::matrix append_column(const GiNaC::matrix& _matrix, const GiNaC::matrix& _col) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r, c + 1);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
        newmat(i, c) = _col(i, 0);
    }
    return newmat;
}

static GiNaC::matrix append_row(const GiNaC::matrix& _matrix, const GiNaC::matrix& _row) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r + 1, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
    for (int j = 0; j < c; j++)
        newmat(r, j) = _row(0, j);
    return newmat;
}

/**
 * Compute a basis of the null space of `mat`. Each column of the returned
 * matrix represents a basis vector.
 */
static GiNaC::matrix nullspace(const GiNaC::matrix& mat) {
    int r = mat.rows(), c = mat.cols();
    
    std::vector<GiNaC::symbol> vars;
    for (int i = 0; i < c; i++)
        vars.push_back(GiNaC::symbol("x" + std::to_string(i)));

    GiNaC::matrix X(c, 1), B(r, 1);
    for (int i = 0; i < c; i++)
        X(i, 0) = vars[i];
    GiNaC::matrix solution = mat.solve(X, B);
    GiNaC::matrix indep_cols(c, 0);
    for (int i = 0; i < c; i++) {
        GiNaC::lst rule;
        for (int j = 0; j < c; j++)
            if (i == j)
                rule.append(vars[j] == 1);
            else
                rule.append(vars[j] == 0);
        
        auto new_indep_cols = append_column(indep_cols, SUBS(solution, rule));
        if (new_indep_cols.rank() > indep_cols.rank())
            indep_cols = new_indep_cols;
    }

    return indep_cols;
}

static int vector_in_difference(const GiNaC::matrix& V, const GiNaC::matrix& W) {
    if (V.rows() == 0)
        return -1;
    
    if (W.rows() == 0)
        return 0;

    int n = W.cols(), found = -1;
    unsigned r = W.rank();

    for (int i = 0; i < V.rows(); i++) {
        if (append_row(W, SUBMATRIX(V, i, 1, 0, n)).rank() > r) {
            found = i;
            break;
        }
    }

    return found;
}


GiNaC::lst jordan::eigenvalues(int* mult) {
    GiNaC::symbol lambda("lambda");
    GiNaC::ex charpoly = GiNaC::factor(A.charpoly(lambda).normal().numer().expand());

    // examine each factor of charpoly
    GiNaC::lst factors, eigenvalues;
    set_factors(charpoly, factors);

    int i = 0;
    for (auto& factor: factors) {
        auto expanded = factor.expand().collect(lambda);
        int degree = expanded.degree(lambda);

        GiNaC::ex coeff1, coeff0, base;
        switch (degree) {
            case 0: // ignore
                break;
            case 1: // a single root
                coeff1 = expanded.coeff(lambda, 1);
                coeff0 = expanded.coeff(lambda, 0);
                eigenvalues.append((-coeff0 / coeff1).normal());
                mult[i++] = 1;
                break;
            default: // non-rational roots, or a multiple root
                if (GiNaC::is_exactly_a<GiNaC::add>(factor))
                    throw std::runtime_error("The matrix has non-rational eigenvalues");
                if (!GiNaC::is_exactly_a<GiNaC::power>(factor))
                    throw std::logic_error("Runs into unexpected case!");
                
                base = factor.begin()->expand().collect(lambda);
                coeff1 = base.coeff(lambda, 1);
                coeff0 = base.coeff(lambda, 0);
                eigenvalues.append((-coeff0 / coeff1).normal());
                mult[i++] = degree;
                break;
        }
    }

    return eigenvalues;
}


GiNaC::matrix jordan::jordan_blocks(const GiNaC::lst& eigvalues, const int* mult, int* num_blocks,
                                    int* block_lambda, int* block_size) {
    int neigs = eigvalues.nops();
    std::vector<int> ranks(N() + 1), diagram(N());
    *num_blocks = 0;
    
    for (int i = 0; i < neigs; i++) {
        int exp_i = mult[i];
        if (exp_i == 1) {
            block_lambda[*num_blocks] = i;
            block_size[*num_blocks] = 1;
            *num_blocks += 1;
        } else {
            GiNaC::matrix B(N(), N()), C;
            for (int j = 0; j < N(); j++)
                for (int k = 0; k < N(); k++)
                    if (j == k)
                        B(j, j) = A(j, j) - eigvalues[i];
                    else
                        B(j, k) = A(j, k);
            
            C = NORMAL(B);
            int ranks_len = 2;
            
            ranks[0] = N();
            ranks[1] = C.rank();
            int j = 0;
            while (ranks[j] > ranks[j + 1] && ranks[j + 1] + exp_i > N()) {
                C = NORMAL(B.mul(C));
                ranks[ranks_len++] = C.rank();
                j++;
            }

            for (j = 0; j < ranks_len - 1; j++)
                diagram[j] = ranks[j] - ranks[j + 1];

            for (j = 1; j <= diagram[0]; j++) {
                int c = 0;
                for (int k = 0; k < ranks_len - 1; k++)
                    c += (diagram[k] >= j);

                block_lambda[*num_blocks] = i;
                block_size[*num_blocks] = c;
                *num_blocks += 1;
            }
        }
    }

    int count = 0;
    for (int i = 0; i < *num_blocks; i++) 
        count += block_size[i];

    if (count != N())
        throw std::logic_error("Sum of block sizes does not agree with size of matrix!");

    GiNaC::matrix J(N(), N());

    count = 0;
    for (int i = 0; i < *num_blocks; i++) {
        for (int j = 0; j < block_size[i]; j++) {
            J(count + j, count + j) = eigvalues[block_lambda[i]];
            if (j < block_size[i] - 1)
                J(count + j, count + j + 1) = 1;
        }
        count += block_size[i];
    }

    return J;
}


GiNaC::matrix jordan::jordan_transform(const GiNaC::lst& eigvalues, int num_blocks, int* block_lambda,
                                       int* block_size) {
    if (N() == 0)
        return GiNaC::matrix(0, 0);
    
    int num_lambda = eigvalues.nops();
    GiNaC::matrix result_mat(N(), N());

    if (num_lambda == N()) {
        GiNaC::matrix B(N(), N());

        for (int i = 0; i < N(); i++) {
            for (int j = 0; j < N(); j++)
                for (int k = 0; k < N(); k++)
                    if (j == k)
                        B(j, j) = A(j, j) - eigvalues[block_lambda[i]];
                    else
                        B(j, k) = A(j, k);

            GiNaC::matrix Y = NORMAL(nullspace(NORMAL(B)));

            if (Y.cols() != 1)
                throw std::logic_error("Runs into unexpected case!");

            for (int j = 0; j < N(); j++)
                result_mat(j, i) = Y(j, 0);
        }

        return result_mat;
    }

    std::vector<int> sizes(N()), counts(N()), written(num_blocks);
    GiNaC::matrix B(N(), N());

    for (int i = 0; i < num_lambda; i++) {
        for (int j = 0; j < N(); j++)
            for (int k = 0; k < N(); k++)
                if (j == k)
                    B(j, j) = A(j, j) - eigvalues[i];
                else
                    B(j, k) = A(j, k);

        B = NORMAL(B);

        int num_sizes = 0;
        for (int j = 0; j < num_blocks; j++) {
            if (block_lambda[j] == i) {
                int k = 0;

                for (k = 0; k < num_sizes; k++) {
                    if (sizes[k] == block_size[j]) {
                        counts[k]++;
                        break;
                    }
                }

                if (k == num_sizes) {
                    sizes[num_sizes] = block_size[j];
                    counts[num_sizes] = 1;
                    num_sizes++;
                }
            }
        }

        GiNaC::matrix Y(N(), N());
        int y_rows = 0;
        for (int j = 0; j < num_sizes; j++) {
            int size = sizes[j];
            int count = counts[j];
            
            if (size == 0)
                throw std::logic_error("Runs into unexpected case!");

            GiNaC::matrix V2 = NORMAL(B.pow(size - 1)),
                          V1 = NORMAL(B.mul(V2)),
                       V1ker = NORMAL(nullspace(V1)).transpose(),
                       V2ker = NORMAL(nullspace(V2)).transpose();

            for (int k = 0; k < count; k++) {
                int rV2ker = V2ker.rows();
                GiNaC::matrix V2kerY(rV2ker + y_rows, N());
                for (int m = 0; m < rV2ker; m++)
                    for (int q = 0; q < N(); q++)
                        V2kerY(m, q) = V2ker(m, q);
                for (int m = 0; m < y_rows; m++)
                    for (int q = 0; q < N(); q++)
                        V2kerY(rV2ker + m, q) = Y(m, q);

                int v_index = vector_in_difference(V1ker, V2kerY);
                if (v_index == -1)
                    throw std::logic_error("Runs into unexpected case!");
                
                int column_offset = 0, output_block = 0;
                while (true) {
                    if (block_lambda[output_block] == i && block_size[output_block] == size && !written[output_block]) {
                        written[output_block] = 1;
                        break;
                    } else {
                        column_offset += block_size[output_block];
                        output_block++;
                    }
                }

                for (int q = 0; q < N(); q++)
                    Y(y_rows, q) = V1ker(v_index, q);

                for (int m = 1; m < size; m++) {
                    auto next_row = NORMAL(B.mul(SUBMATRIX(Y, y_rows + m - 1, 1, 0, N()).transpose()));
                    for (int q = 0; q < N(); q++)
                        Y(y_rows + m, q) = next_row(q, 0);
                }
                y_rows += size;

                for (int m = 0; m < size; m++)
                    for (int l = 0; l < N(); l++)
                        result_mat(l, column_offset + m) = Y(y_rows - 1 - m, l);
            }
        }
    }

    return result_mat;
}


