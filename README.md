# Matrix Decomposition Methods (Numerical Linear Algebra)
Implementation and analysis of core numerical matrix decomposition algorithms, with a focus on understanding how linear algebra theory translates into practical computation and stability considerations relevant to data science and machine learning.

**Overview:**
- Implemented LU, QR (Classical Gram–Schmidt), SVD, and Schur decompositions from scratch in MATLAB
- Built the Schur decomposition using an iterative QR algorithm, reusing the custom QR implementation
- Constructed SVD via the custom Schur decomposition applied to A^T*A (not built-in functions)
- Evaluated reconstruction error, numerical stability, and runtime behavior across increasing matrix sizes
- Examined the impact of conditioning and floating-point precision on algorithm performance
- Presented results and interpretations in a structured technical report, discussing when and why the custom functions behaved poorly.

**Tools:**
- MATLAB
- Advanced/Numerical Linear Algebra
- Error/Sensitivity analysis

**Files:**

src/ — MATLAB implementations of LU, QR, Schur, and SVD

experiments/ — numerical error, stability, and runtime analyses

report/ — technical report (PDF)

