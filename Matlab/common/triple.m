% Function implementing the scalar triple product:
% Statistical Optimizatopn for Geometric Computation, Kanatani p.32

function product = triple(x1,x2,x3)

product = dot(cross(x1,x2), x3);