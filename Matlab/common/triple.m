%% Function implementing the scalar triple product Kanatani p.32

function product = triple(x1,x2,x3)

product = dot(cross(x1,x2), x3);