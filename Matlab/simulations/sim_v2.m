clearvars
clc

%% Initialize
depth_range = 10;
N_points = 10;

C1_P = zeros(3, N_points); %% Euclidean (euclidean values)
C1_pix = zeros(2, N_points); %% Pixel (pixel values)
C1_h = zeros(3, N_points); %% Homogeneous (normalized euclidean values)
C1_b = zeros(3, N_points); %% Bearings (normalized euclidean values)
C1_d = zeros(1, N_points); %% Depths

C2_P = zeros(3, N_points); %% Euclidean (euclidean values)
C2_pix = zeros(2, N_points); %% Pixel (pixel values)
C2_h = zeros(3, N_points); %% Homogeneous (normalized euclidean values)
C2_b = zeros(3, N_points); %% Bearings (normalized euclidean values)

%% Intrinsics
camera.fc = [400;400];
camera.cc = [320;240];
camera.dx = 640;
camera.dy = 480;

%% Create random motion parameters
% Rotations
% reference R1
C1_R_C1 = eye(3);
% reference R2
C2_R_C1 = generate_random_rotation();

% Translations
% reference t1
C1_P_C1 = [0; 0; 0];
% reference t2
C2_P_C1 = [6 * rand; 2 * rand ; 1 * rand];

%% Populate the variables
% Pay attention to the fact that there will always be an error when
% transforming between frames of reference. This is why the 3D point as
% seen in camera 1 (C1_p) will not match exactly to the same 3D point after 
% transforming to C2 (C2_p) and back (new_C1_p):
% C1_p -> C2_p = (R * C1_p + t) -> new_C1_p ~= (R'*C2_p - R'*t)
% C1_p does not match exactly to the new_C1_p!
noise_gain = 0.1;
for i = 1:N_points
    random_x_pixel = round(rand * camera.dx);
    random_y_pixel = round(rand * camera.dy);
    random_depth = 5 + round(rand * depth_range);
        
    pixel = [random_x_pixel; random_y_pixel];
    homogeneous = [(pixel - camera.cc) ./ camera.fc ; 1];
    bearing = homogeneous / norm(homogeneous);
    euclidean = bearing * random_depth;
    
    % Noise in measurements 1
    g_noise_1 = [noise_gain*rand; noise_gain*rand; 0];
    
    homogeneous_noisy = homogeneous + g_noise_1;
    bearing_noisy = homogeneous_noisy / norm(homogeneous_noisy);
    
    C1_P(:,i) = euclidean;
    C1_d(:,i) = random_depth;
    C1_h(:,i) = homogeneous_noisy;
    C1_b(:,i) = bearing_noisy;
    C1_pix(:,i) = pixel;
            
%     g_noise_2 = noise_gain*[rand; rand; rand];
    
    C2_euclidean = C2_R_C1' * (euclidean) - C2_R_C1 * C2_P_C1;
    C2_homogeneous = C2_euclidean / C2_euclidean(3);
    C2_bearing = C2_homogeneous / norm(C2_homogeneous);
    C2_pixel = C2_homogeneous(1:2) .* camera.fc + camera.cc;
    
    C2_P(:,i) = C2_euclidean;
    C2_h(:,i) = C2_homogeneous;
    C2_b(:,i) = C2_bearing;
    C2_pix(:,i) = round(C2_pixel);
end

%% Essential from known motion parameters
E_bearings =  skew_sym(C2_P_C1) * C2_R_C1;
epipolar_b = diag( C2_b' * E_bearings * C1_b );
% % Test epipolar constraint. This is a test for the validity of the 
% % initialization operations and should never assert when noise_gain = 0!!!
% for i = 1:N_points
%     assert(norm( C2_b(:,i)' * E_bearings * C1_b(:,i)) < 1e-10)
% end

%% OpenGV Kneip
% Run the 5 point algorithm from Nister to get motion parameters. Remember
% the function requires bearing parameters!
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );

%% Find Rotation and Translation {R,t} for the second camera---------------
R2 = R_t(:,1:3);
t2 = R_t(:,4);

%% Estimate Essential matrix
E_estimated = skew_sym(t2) * R2; % E = skew_sym(t) * R
epipolar = diag( C2_b' * E_estimated * C1_b );
% Test epipolar constraint. This may assert for very low threshold, it is
% an estimation with noisy data after all.
% for i = 1:N_points
% %     assert(norm( C2_h(:,i)' * E_estimated * C1_h(:,i)) < 1e-1)
%     epipolar(i) = norm( C2_h(:,i)' * E_estimated * C1_h(:,i) );
% end

%% Kanatani correction
reps = 1;
[v_k, u_k, Dv_k, Du_k, keepers, new_R2, new_t2] = correct_points_kanatani_iterative(C1_b, C2_b, R2, t2, E_estimated, reps);

% Test epipolar constraint with Kanatani correction.
new_Essential = skew_sym(t2) * R2;
epipolar_k = diag( u_k' * new_Essential * v_k );

disp(['Epipolar using Kanatani points: ' num2str(norm(epipolar_k))])
disp(['Epipolar using E_estimated:     ' num2str(norm(epipolar))])
disp(['Epipolar using E_bearings:      ' num2str(norm(epipolar_b))])

% %% Intrinsics matrix
% K = [camera.fc(1)      0       camera.cc(1);
%           0       camera.fc(2) camera.cc(2);
%           0            0            1      ];
% 
% F_bearings = (inv(K)') * E_bearings * inv(K);
% 
% %% Epipolar constraint for Fundamental matrix
% for i = 1:N_points
%     assert(norm( [C2_pix(:,i);1]' * F_bearings * [C1_pix(:,i);1]) < 1e-10)
% end

%% Triangulate kanatani
R1 = eye(3);
t1 = [0 0 0]';

for i=1:size(v_k,2)
    % No Kanatani correction reconstruction
    Z1(i) = dot(cross(t2,R2*C2_h(:,i)), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
    Z2(i) = dot(cross(C1_h(:,i),t2), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
    % Kanatani correction reconstruction
%     Zk1(i) = dot(cross(t2,R2*u_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
%         norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
%     
%     Zk2(i) = dot(cross(t2,v_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
%         norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
%            
%     n = cross(t2,v_k(:,i))/norm(cross(t2,v_k(:,i)));
%     m = cross(n,R2*u_k(:,i));
%     
%     DZk(i) = -(dot(m, Zk1(i)*Dv_k(:,i) - Zk2(i)*R2*Du_k(:,i))) / dot(m,v_k(:,i));
    
    % Kanatani with updated motion aprameters
    Zk1(i) = dot(cross(t2,R2*u_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
    
    Zk2(i) = dot(cross(t2,v_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
           
    n = cross(t2,v_k(:,i))/norm(cross(t2,v_k(:,i)));
    m = cross(n,R2*u_k(:,i));
    
    DZk(i) = -(dot(m, Zk1(i)*Dv_k(:,i) - Zk2(i)*R2*Du_k(:,i))) / dot(m,v_k(:,i));
    
    r1(1:3,i) = ( Zk1(i) + DZk(i) )*v_k(:,i);
    r2(1:3,i) = Z1(i)*C1_h(:,i);
    
    r1_n(1:3,i) = r1(1:3,i) / norm( r1(1:3,i) );
    r2_n(1:3,i) = r2(1:3,i) / norm( r2(1:3,i) );
end

%% Reprojection error
error_k = reprojection_error(R1, t1, R2, t2, r1, C1_b(:,keepers), C2_b(:,keepers));
error   = reprojection_error(R1, t1, R2, t2, r2, C1_b(:,keepers), C2_b(:,keepers));
disp(['Kanatani reprojection error: ' num2str(error_k)])
disp(['Plain reprojection error:    ' num2str(error)])

%% Rotation error (paper Metrics for 3D Rotations: Comparison and Analysis)
rot_error_k = norm(eye(3) - (C2_R_C1*R2'),'fro');
rot_error = norm(eye(3) - (C2_R_C1*R2'),'fro');
disp(['Kanatani rotation error: ' num2str(rot_error_k)])
disp(['Plain rotation error:    ' num2str(rot_error)])

%% Translation error
trans_error_k = acos( (C2_P_C1' * new_t2) / (norm(C2_P_C1) * norm(new_t2)) );
trans_error = acos( (C2_P_C1' * t2) / (norm(C2_P_C1) * norm(t2)) );
disp(['Kanatani translation error: ' num2str(trans_error_k)])
disp(['Plain translation error:    ' num2str(trans_error)])

%% Express all points w.r.t. their camera for visualization purposes
vis_pts1 = []; vis_pts2 = [];
for i=1:size(C1_h,2)
    vis_pts1(1:3,i) = C1_h(:,i);
    vis_pts2(1:3,i) = R2*C2_h(:,i) + t2;
    vis_r1(1:3,i) = r1(1:3,i);
    vis_r2(1:3,i) = r2(1:3,i);
end

figure, xlabel('x'), ylabel('y'), zlabel('z')
axis equal, grid on, hold on
plot3(t1(1), t1(2), t1(3), 'rx')
plot3(t2(1), t2(2), t2(3), 'bx')
for i=1:size(C1_h,2)
    plot3([t1(1) vis_pts1(1,i)], [t1(2) vis_pts1(2,i)], [t1(3) vis_pts1(3,i)], 'b')
    plot3([t2(1) vis_pts2(1,i)], [t2(2) vis_pts2(2,i)], [t2(3) vis_pts2(3,i)], 'g')
    plot3(C1_P(1,i), C1_P(2,i), C1_P(3,i), 'b^')
    plot3(vis_r1(1,i), vis_r1(2,i), vis_r1(3,i), 'b.')
    plot3(vis_r2(1,i), vis_r2(2,i), vis_r2(3,i), 'r.')
end
visualize_camera_pose(camera, R1, t1)
visualize_camera_pose(camera, R2, t2)