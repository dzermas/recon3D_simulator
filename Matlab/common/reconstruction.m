function [r1_filtered, r2_filtered, motion_params] = reconstruction(C1_h, C2_h)
%% Variables
% % IN:
% % C*_homogeneous == 2D in 3Dimensions, not normalized aka norm != 1
% % OUT:
% % r1_filtered : 3D points with Kanatani correction
% % r2_filtered : 3D points without Kanatani correction
% % motion_params : structure that holds all motion parameters (R1, h1, R2, h2, Essential matrix)

% C*_bearing == normalized 3 homogeneous vector
for i=1:size(C1_h,2)
    C1_b(1:3,i) = C1_h(:,i) / norm(C1_h(:,i));
    C2_b(1:3,i) = C2_h(:,i) / norm(C2_h(:,i));
end

%% Try OpenGV Kneip to find motion parameters
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );

%% Find Rotation and Translation {R,t} for the second camera---------------
R2 = R_t(:,1:3);
h2 = R_t(:,4);
Essential = skew_sym(h2) * R2; % E = skew_sym(t) * R

% 1-2 reps are enough, the error doesn't drop any more
for reps=1:2
    
    %% Correct points to accurately validate the epipolar constraint-----------
    % [v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani(C1_h, C2_h, Essential);
    
    [v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani_iterative(C1_h, C2_h, R2, h2, Essential, reps);
    
    %% Triangulate kanatani
    R1 = eye(3);
    h1 = [0 0 0]';
    
    Z.Z1 = []; Z.Zk1 = []; Z.Zk2 = []; % Z.Z2 = [];
    for i=1:size(v_k,2)
        % No Kanatani correction reconstruction
        Z1 = dot(cross(h2,R2*C2_h(:,i)), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
            norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
        
        Z2 = dot(cross(C1_h(:,i),h2), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
            norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
        
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
        
        Z.Z1(i) = Z1;
        %             Z.Z2(counter) = Z2;
        Z.Zk1(i) = Zk1;
        Z.Zk2(i) = Zk2;
    end

    %% Remove reconstructed points at infinity or behind the camera
    counter_infinity = 1;
    for i=1:length(Z.Z1)
        if Z.Z1(i) < 2*mean(Z.Z1) && Z.Z1(i) > 0
            r1_filtered(1:3,counter_infinity) = r1(:,i);
            r2_filtered(1:3,counter_infinity) = r2(:,i);
            counter_infinity = counter_infinity + 1;
            pass(i) = true;
        else
            pass(i) = false;
        end
    end

    %% Reprojection error
    error_k = reprojection_error(R1, h1, R2, h2, r1_filtered, C1_b(:,(keepers & pass)), C2_b(:,(keepers & pass)));
    error   = reprojection_error(R1, h1, R2, h2, r2_filtered, C1_b(:,pass), C2_b(:,pass));
    disp(['Kanatani reprojection error: ' num2str(error_k)])
    disp(['Plain reprojection error:    ' num2str(error)])
end

motion_params.R1 = R1;
motion_params.h1 = h1;
motion_params.R2 = R2;
motion_params.h2 = h2;
motion_params.E  = Essential;