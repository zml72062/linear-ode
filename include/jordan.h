#ifndef JORDAN_H
#define JORDAN_H


/**
 * jordan.h - Symbolic Jordan decomposition with GiNaC.
 */

#include <ginac/ginac.h>


struct jordan {
    /**
     * Compute symbolically the Jordan canonical form of `mat`.
     */
    jordan(const GiNaC::matrix& mat): A(mat) {
        if (mat.rows() != mat.cols())
            throw std::runtime_error("Not a square matrix");
    }

    GiNaC::matrix A;

    int N() { return A.rows(); }

    /**
     * Return a list of eigenvalues, and put their respective multiplicities
     * in the array `mult`. The caller is responsible for allocating memory 
     * for `mult`. 
     * 
     * When the function fails to express all eigenvalues in rational form,
     * it will raise a `std::runtime_error`.
     */
    GiNaC::lst eigenvalues(int* mult);

    /**
     * Build Jordan canonical form given the list of eigenvalues `eigvalues`
     * and their respective multiplicities `mult`.
     * 
     * The function sets `*num_blocks` to the number of Jordan blocks, 
     * `block_lambda` to a sequence whose i-th entry represents the index into
     * `eigvalues` of eigenvalue corresponding to the i-th Jordan block, and
     * `block_size` to a sequence saving the block sizes. Finally, it returns
     * the block-diagonal Jordan canonical matrix.
     */
    GiNaC::matrix jordan_blocks(const GiNaC::lst& eigvalues, const int* mult, int* num_blocks,
                                int* block_lambda, int* block_size);

    /**
     * Give matrix P such that P^(-1)AP is a Jordan canonical matrix. The
     * caller should provide precomputed list of eigenvalues `eigvalues`,
     * and information about the Jordan canonical matrix, namely `num_blocks`,
     * `block_lambda` and `block_size`. Please check `jordan_blocks()` for
     * the definition of those fields.
     */
    GiNaC::matrix jordan_transform(const GiNaC::lst& eigvalues, int num_blocks, int* block_lambda,
                                   int* block_size);

};


#endif // JORDAN_H

