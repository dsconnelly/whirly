# whirly

`whirly` provides easy-to-use Python classes for the spectral solutions of two-dimensional partial differential equations, especially those related to fluid dynamics. For example, `whirly` includes a solver for the two-dimensional vorticity equation
$$\frac{\partial \zeta}{\partial t} + \mathbf{u} \cdot \nabla \zeta = \frac{1}{\text{Re}}\nabla^2 \zeta$$
where $\text{Re}$ is the Reynolds number. For time-stepping, `whirly` implements a fourth-order integrating factor method as outlined in [Yang et al. (2021)](https://www.sciencedirect.com/science/article/pii/S002199912030766X).

See `example.ipynb` for a quick demonstration of `whirly` in action.
