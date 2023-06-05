# CSE306_assignment2
A free-surface 2D fluid solver using incompressible Euler's equations.

This project includes the following :
- A Voronoï diagram using Voronoï Parallel Linear Enumeration (Sec. 4.3.3, lab 6) with Sutherland-Hodgman polygon clipping algorithm (Sec. 4.2, lab 6).
- Extending this to a Power diagram (Sec. 4.4.3, lab 7).
- Optimizing the weights of the power diagram using LBFGS (Sec 4.4.4, lab 7)
- Implementing de Gallouet-Mérigot incompressible Euler scheme that prescribes each fluid cell area and the sum of air cell areas to given values, and add a spring force from each fluid particle to their Laguerre's cell centroid (Sec. 5.4, lab 8).
