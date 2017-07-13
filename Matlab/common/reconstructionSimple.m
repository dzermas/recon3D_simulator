function [r1, r2] = reconstructionSimple(C1_h, C2_h, R2, h2)
%% Variables
% % IN:
% % C*_homogeneous == 2D in 3Dimensions, not normalized aka norm != 1
% % R2,h2 : transformation parameters for second camera
% % OUT:
% % r1_filtered : 3D points with Kanatani correction
% % r2_filtered : 3D points without Kanatani correction



%% Correct points to accurately validate the epipolar constraint-----------
Essential = skew_sym(h2) * R2; 
[v_k, u_k, Dv_k, Du_k] = correct_points_kanatani(C1_h, C2_h, Essential);

%% Triangulate kanatani
Z.Z1 = []; Z.Zk1 = []; Z.Zk2 = [];
for i=1:size(v_k,2)
    % No Kanatani correction reconstruction
    Z1 = dot(cross(h2,R2*C2_h(:,i)), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
%     Z2 = dot(cross(C1_h(:,i),h2), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
%         norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
    % Kanatani correction reconstruction
    Zk1 = dot(cross(h2,R2*u_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
    
    Zk2 = dot(cross(h2,v_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
    
    n = cross(h2,v_k(:,i))/norm(cross(h2,v_k(:,i)));
    m = cross(n,R2*u_k(:,i));
    
    DZk = -(dot(m, Zk1*Dv_k(:,i) - Zk2*R2*Du_k(:,i))) / dot(m,v_k(:,i));
    
    r1(1:3,i) = (Zk1 + DZk)*v_k(:,i);
    r2(1:3,i) = Z1*C1_h(:,i);
end