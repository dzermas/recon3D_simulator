function [v_k, u_k, Dv_k, Du_k, keepers, R2_k, t2_k] = correct_points_kanatani_iterative(C1_b, C2_b, R2, t2, E_estimated, reps)

for i=1:reps
    [v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani(C1_b, C2_b, E_estimated);
    
    epipolar_k = diag( C2_b' * E_estimated * C1_b );
    
    %% Need to implement {DR | Dt} update here!!!
    R2_k = rotation_correction_katnatani(C1_b, C2_b, R2, t2, epipolar_k);
    
    t2_k = translation_correction_katnatani(C1_b, C2_b, R2, t2, epipolar_k);
    
%     for j=1:size(v_k,2)
%         v_k_b(1:3,j) = v_k(:,j) / norm(v_k(:,j));
%         u_k_b(1:3,j) = u_k(:,j) / norm(u_k(:,j));
%     end
%     
%     E_corrected = skew_sym(t2_k) * R2;
end
sum(keepers);