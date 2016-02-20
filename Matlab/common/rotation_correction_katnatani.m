function R2_c = rotation_correction_katnatani(v_k, u_k, R2, t2, epipolar_k)

b_v = [];
for i=1:size(v_k,2)
    b = dot( v_k(:,i), R2*u_k(:,i) ) * t2 - dot( t2, R2*u_k(:,i) ) * v_k(:,i);
    b_v = [b_v b];
end

Domega = - inv(b_v*b_v') * b_v * epipolar_k;

DR = skew_sym(Domega) * R2;

R2_c = R2 + DR;