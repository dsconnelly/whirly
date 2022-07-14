whirly
======

whirly provides easy-to-use Python classes for the spectral solutions of two-dimensional partial differential equations, especially those related to fluid dynamics. For example, whirly includes a solver for the two-dimensional vorticity equation

.. math:: \frac{\partial \zeta}{\partial t} + \mathbf{u} \cdot \nabla \zeta = \frac{1}{\textrm{Re}}\nabla^2 \zeta

where :math:`\textrm{Re}` is the Reynolds number. For time stepping, whirly implements a fourth-order integrating factor method as outlined in `Yang et al. (2021) <https://www.sciencedirect.com/science/article/pii/S002199912030766X>`_.

See the `example notebook <https://github.com/dsconnelly/whirly/blob/main/example.ipynb>`_ for a quick demonstration of whirly in action.
