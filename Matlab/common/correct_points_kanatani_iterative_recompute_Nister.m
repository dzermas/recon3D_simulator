function [v_k, u_k, Dv_k, Du_k, keepers, R2_k, h2_k] = ...
    correct_points_kanatani_iterative_recompute_Nister(v_k, u_k, Essential, reps)

for i=1:reps
    [v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani(v_k, u_k, Essential);
    
    for j=1:size(v_k,2)
        v_k_b(1:3,j) = v_k(:,j) / norm(v_k(:,j));
        u_k_b(1:3,j) = u_k(:,j) / norm(u_k(:,j));
    end
    
    R_t_k = opengv ( 'fivept_nister_ransac', v_k_b, u_k_b );
    
    R2_k = R_t_k(:,1:3);
    h2_k = R_t_k(:,4);
    Essential = skew_sym(h2_k) * R2_k;
end
sum(keepers);