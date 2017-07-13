function P = singlePointReconstruction(C1_h, C2_h, R, t)


Z1 = dot(cross(t,R*C2_h), cross(C1_h,R*C2_h)) / ...
        norm(cross(C1_h,R*C2_h))^2;
    
P = Z1*C1_h;