# whirly

whirly provides easy-to-use Python classes for the spectral solution of two-dimensional partial differential equations, especially those related to fluid dynamics. The core of the package is the abstract `PseudospectralSolver` class, which is easily extensible to a variety of PDEs with linear and nonlinear terms. For time stepping, whirly implements a fourth-order integrating factor method as outlined in [Yang et al. (2021)](https://www.sciencedirect.com/science/article/pii/S002199912030766X), accessible through the `IFRK4Integrator` class.

See [`example.ipynb`](https://github.com/dsconnelly/whirly/blob/main/example.ipynb) for a quick demonstration of whirly used to solve the two-dimensional [vorticity equation](https://en.wikipedia.org/wiki/Vorticity_equation).
