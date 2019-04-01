# CDFSampler

Sampling discrete 1D and 2D CDFs "the right way".

A pdf is a function, which we are here given discretised in an array of non-negative values f_1, ..., f_N of sum equal to 1. Or in 2D, a table of double-indexed values.

Optionally, we are also given an array of positions x_1, ..., x_N; otherwise only the edge positions and N. For simplicity, we will omit the 2D case from now on, but it is implemented as well.

We can view f_1, ..., f_N as a distribution of Dirac masses at x_1, ..., x_N. For the purpose of sampling, we need to low-pass filter these masses back into a function.

The simplest of such filters is the box filter. However, the width of the filtering box is allowed to vary (as long as the filter's L^1 norm remains 1). We chose for the support of these filters the Voronoi region around each data point x_1,...,x_N (prolonged by symmetry across the edges).

In other words, the filtered function is piecewise constant on [x_1, (x_1+x_2)/2], [(x_1+x_2)/2,(x_2+x_3)/2], ..., [(x_{N-1}+x_N)/2, x_N], such that the integral on the i-th interval is f_i.

The cdf, defined as F(x) = \int_x_1^x f(x) dx, is thus piecewise affine on the same intervals. Since we know that F(x_1) = 0, we only store values at the N positions (x_1+x_2)/2, ..., (x_{N-1}+x_N)/2, x_N. There are N values to store for F, same array size as f which is convenient but remember that the positions at which these functions were sampled are different. If we work with a non-uniform distributions of positions we also store the "staggered" positions for F.

Once this is understood the code should be self-explanatory. See the Doxygen'd documentation.

class CDFSampler holds dimension-agnostic data.
Dimensions 1 and 2 are supported in CDF1Sampler and CDF2Sampler respectively, each declined in
CDFxSamplerUniform and CDFxSamplerNonUniform.
