function d = distance_between_lines(p1,p2,n1,n2)

if cross(p1,p2) == 0
    d = norm(n1/norm(p1) - n2/norm(p2));
else
    d = (dot(p1,n2) + dot(p2,n1)) / norm(cross(p1,p2));
end