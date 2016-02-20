function t2_c = translation_correction_katnatani(v_k, u_k, R2, t2, epipolar_k)

a_v = [];
for i=1:size(v_k,2)
    a = cross( v_k(:,i), R2*u_k(:,i) );
    a_v = [a_v a];
end

Dt = - inv(a_v*a_v') * a_v * epipolar_k;

t2_c = t2 + Dt;