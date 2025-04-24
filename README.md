# Solving Homogeneous Linear Differential Systems with Power Series

This repository implements an algorithm presented by Moulay Barkatou ([link](https://www.impan.pl/~slawek/pisa/Barkatou_p.pdf)) to solve a homogeneous linear differential system around either an analytic point or a regular singularity.

## Dependencies

The program uses `ginac` as its frontend, while exploits either `flint` or `ginac` as its backend. One must have libraries `flint` and `ginac` installed on the computer.

## Usage

Run `make` to build the program. Please check `example.cpp` and `symexample.cpp` to see example usage. 

## Introduction

A homogeneous linear differential system can be written as $\frac{dY}{dx} = A(x)Y(x)$, where $A(x)\in\mathbb{R}^{N\times N}$ and $Y(x)\in\mathbb{R}^{N}$. If $A(x)$ is analytic around $x=x_0$, the point $x_0$ is called an **analytic point** of the equation. If $A(x)$ is singular at $x=x_0$, but an arbitrary solution $Y(x)$ grows no faster than a polynomial of $|x-x_0|^{-1}$ as $x$ approaches $x_0$, then the point $x_0$ is called a **regular singularity**. In other cases, $x_0$ is called an **irregular singularity**.

It is proved that an analytic point of the equation must be an analytic point of the solution. Therefore, around an analytic point $x_0$ we can Taylor expand $Y(x)$ as
$$Y(x)=\sum_{k=0}^\infty Y_k(x-x_0)^k$$
If we analogously expand $A(x)$ into a Taylor series, we can determine $Y_k$ iteratively, using the formula
$$Y_{k+1}=\frac{1}{k+1}\left(\sum_{j=0}^kA_{k-j}Y_j\right)$$
Now $Y_0$ is the arbitrary constant to be determined by initial condition at $x=x_0$.

If $x_0$ is a regular singularity, one may show that a complete set of linearly independent solutions can be found by considering solutions of the form
$$Y(x)=x^\lambda \sum_{p=0}^{t-1}\frac{\log^px}{p!}\sum_{k=0}^\infty Y_k^{(p)}x^k$$
where $\lambda\in\mathbb{C}$ is a complex constant and $t$ is a finite positive integer. 

However, determining $Y_k^{(p)}$ directly is usually hard. An effective process to tackle regular singularities is stated by Moulay Barkatou ([link](https://www.impan.pl/~slawek/pisa/Barkatou_p.pdf)). It is divided into three stages. The first stage is to apply a reduction procedure (known to be explored by Moser) to reduce the pole order of $A(x)$ at $x=x_0$ to one. After the reduction one may write
$$A(x)=\frac{1}{x}\sum_{k=0}^\infty A_kx^k$$
and the recurrence relation for $Y_k^{(p)}$ goes as
$$[(\lambda+k)I-A_0]Y_k^{(p)}+Y_k^{(p+1)}1[p\ne t-1]=\sum_{l=0}^{k-1}A_{k-l}Y_l^{(p)}$$
When $k=0$, the above equation reduces to
$$(\lambda I-A_0)Y_0^{(t-1)}=0$$
and $$(\lambda I-A_0)Y_0^{(p)}=-Y_0^{(p+1)}$$
for each $p=0,1,\ldots, t-2$. Now it is clear that the allowed values for $\lambda$ must be eigenvalues of $A_0$, and the leading vectors $Y_0^{(p)}$ are decided by the generalized eigenvectors (a.k.a. Jordan chains) corresponding to $\lambda$.

When one tries to use the recurrence relation to determine $Y_k^{(p)}$ for $k>0$, another problem may occur: the matrix $(\lambda+k)I-A_0$ can be singular for a nonzero $k$, if $\lambda+k$ is also an eigenvalue of $A_0$. To deal with this problematic case, we need to add an additional stage before solving the recurrence relation, in which we reduce $A(x)$ further to make $A_0$ have no eigenvalues that differ from one another by nonzero integers. After that, the coeffcient vectors can be safely determined.

## About the program

From the introduction, one already sees that computing Jordan decomposition is a necessary step when trying to solve an equation around its regular singularities. However, Jordan decomposition is known to be numerically unstable. Therefore, our program prefers performing Jordan decompositions symbolically. 

Actually, we provide two implementations of the solver. The main difference between the two implementations is in how they do Jordan decompositions. One version uses library functions provided by `flint` to do Jordan decomposition. Those library functions are faster, and are guaranteed to give accurate symbolic results when the leading coefficient matrix $A_0$ has rational eigenvalues only. In this case, this implementation will output series solutions with all-rational coefficients. Otherwise, when at least one eigenvalue of $A_0$ is irrational, the corresponding solution will have numerical coefficients. In the latter case, the solution is prone to numerical errors, and the program may even fail to give a solution (when there are highly degenerate roots from a high-degree characteristic polynomial).

Another implementation is purely based on `ginac` to do Jordan decomposition. This has the advantage of always producing accurate symbolic results, and being able to generalize to cases in which the differential equation relies on an additional parameter. However, since it is generally impossible to solve polynomial equations symbolically, this implementation does not guarantee to return a result successfully. Specifically speaking, it will fail and throw an error when at least one eigenvalue of $A_0$ cannot be expressed as a rational function with rational coefficients. Furthermore, computation would be slower compared with the `flint` version.


