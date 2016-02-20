%% 3D line reconstruction simulator
clearvars
clc

%% Initialize
depth_range = 5;
% Need at least 8 points - 4 lines for an estimation of the motion
% parameters
N_lines = 4;
N_points = N_lines * 2;

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
% Generate 3D points
% You can select the gain of Gaussian noise here!
noise_gain = 0.01;
for i = 1:N_points
    % Random pixels on image plain 1...
    random_x_pixel = round(rand * camera.dx);
    random_y_pixel = round(rand * camera.dy);
    random_depth = 5 + round(rand * depth_range);
    
    pixel = [random_x_pixel; random_y_pixel];
    % bring them in homogeneous coordinates...
    homogeneous = [(pixel - camera.cc) ./ camera.fc ; 1];
    % and then in bearing...
    bearing = homogeneous / norm(homogeneous);
    % which then becomes a 3D point with random depth!
    euclidean = bearing * random_depth;
    
    % Save camera 1 values!
    C1_P(:,i) = euclidean;
    C1_d(:,i) = random_depth;
    C1_h(:,i) = homogeneous;
    C1_b(:,i) = bearing;
    C1_pix(:,i) = pixel;
    
    % Transform the 3D points and prepare projecting them on the second
    % camera...
    g_noise = noise_gain*[rand; rand; rand];
    C2_euclidean = C2_R_C1' * (euclidean + g_noise) - C2_R_C1' * C2_P_C1;
    % project them by making them homogeneous...
    C2_homogeneous = C2_euclidean / C2_euclidean(3);
    % and then compute bearings...
    C2_bearing = C2_homogeneous / norm(C2_homogeneous);
    % and pixels just in case we need them!
    C2_pixel = C2_homogeneous(1:2) .* camera.fc + camera.cc;
    
    % Save camera 2 values!
    C2_P(:,i) = C2_euclidean;
    C2_h(:,i) = C2_homogeneous;
    C2_b(:,i) = C2_bearing;
    C2_pix(:,i) = C2_pixel;
end

%% Generate 3D lines from the 3D points and associate the corresponding 2D lines
for i=1:N_lines
    % Create 3D lines by connecting randomly the generated 3D points
    ind = datasample([1:size(C1_P,2)], 2, 'Replace', false);
    line_3D(1:6,i) = [C1_P(:,ind(1)); C1_P(:,ind(2))];
    % Pair up the lines between the two image plains
    line_2D_1(1:6,i) = [C1_h(:,ind(1)); C1_h(:,ind(2))];
    line_2D_2(1:6,i) = [C2_h(:,ind(1)); C2_h(:,ind(2))];
end

%% Estimate Motion Parameters
R1 = eye(3);
t1 = [0;0;0];
% Use the OpenGV library 5 point Niester algorithm
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );
R2 = R_t(:,1:3);
t2 = R_t(:,4);
% R2 = C2_R_C1;
% t2 = C2_P_C1;

%% Reconstruct 3D lines
p = zeros(3, N_lines);
m = zeros(3, N_lines);
rh = zeros(3, N_lines);
n1 = zeros(3, N_lines);
n2 = zeros(3, N_lines);
L = zeros(3, N_lines);
for i=1:N_lines
    init_n1 = cross(line_2D_1(1:3,i), line_2D_1(4:6,i));
    init_n2 = cross(line_2D_2(1:3,i), line_2D_2(4:6,i));
    % Normalized to norm(n,2) = 1
    % normals expressed in the world coordinates (in our case, both 
    % cameras at t = [0 0 0]' and R = eye(3))
    n1(:,i) = init_n1/norm(init_n1,2);
    n2(:,i) = init_n2/norm(init_n2,2);

    % {m,rh}-representation of the 3D line (see Kanatani 4.2.2 and 6.4.1)
    p(:,i) = cross(n1(:,i), R2*n2(:,i)) / dot(t2, R2*n2(:,i));
    m(:,i) = p(:,i) / norm(p(:,i));
    rh(:,i) = cross(p(:,i), n1(:,i)) / norm(p(:,i))^2;
    
    % Begin (L1) and end (L2) of the reconstructed line should be
    % determined by the two points that construct the line. Under noise,
    % the parameters of the 3D line will help correct the 3D position of
    % the end-points.
    Z1(i) = dot(cross(C2_P_C1,R2*line_2D_2(1:3,i)), cross(line_2D_1(1:3,i),R2*line_2D_2(1:3,i))) / ...
        norm(cross(line_2D_1(1:3,i),R2*line_2D_2(1:3,i)))^2;
    
    Z2(i) = dot(cross(C2_P_C1,R2*line_2D_2(4:6,i)), cross(line_2D_1(4:6,i),R2*line_2D_2(4:6,i))) / ...
        norm(cross(line_2D_1(4:6,i),R2*line_2D_2(4:6,i)))^2;
    
    % 3D points
    P1(:,i) = Z1(i)*line_2D_1(1:3,i);
    P2(:,i) = Z2(i)*line_2D_1(4:6,i);
    % The line parameters can be used as constraints to refine 3D points
    L1(:,i) = rh(:,i) - m(:,i);
    L2(:,i) = rh(:,i) + m(:,i);
                
    % Project the 3D points (P) on the 3D line (A,B)
    % A + dot(AP,AB) / dot(AB,AB) * AB
    Pr1(:,i) = L1(:,i) + dot(P1(:,i)-L1(:,i), L2(:,i)-L1(:,i)) ...
        / dot(L2(:,i)-L1(:,i),L2(:,i)-L1(:,i)) * (L2(:,i)-L1(:,i));
    Pr2(:,i) = L1(:,i) + dot(P2(:,i)-L1(:,i), L2(:,i)-L1(:,i)) ...
        / dot(L2(:,i)-L1(:,i),L2(:,i)-L1(:,i)) * (L2(:,i)-L1(:,i));
end

%% Express all points w.r.t. their camera for visualization purposes
% Initialize points to visualize
vis_line_2D_1 = []; vis_line_2D_2 = [];
% Select whether to visualize camera positions in normalized (C*_R_C1 and 
% C*_P_C1) or real (R* and t*) coordinates
vis_R1 = R1;
vis_R2 = R2;
vis_t1 = t1;
vis_t2 = t2;
for i=1:N_lines
    vis_line_2D_1(1:3,i) = vis_R1*line_2D_1(1:3,i) + vis_t1;
    vis_line_2D_1(4:6,i) = vis_R1*line_2D_1(4:6,i) + vis_t1;
    vis_line_2D_2(1:3,i) = vis_R2*line_2D_2(1:3,i) + vis_t2;
    vis_line_2D_2(4:6,i) = vis_R2*line_2D_2(4:6,i) + vis_t2;
    vis_n1(1:3,i) = vis_R1*n1(:,i) + vis_t1;
    vis_n2(1:3,i) = vis_R2*n2(:,i) + vis_t2;
end

figure, xlabel('x'), ylabel('y'), zlabel('z')
axis equal, grid on, hold on
plot3(vis_t1(1), vis_t1(2), vis_t1(3), 'rx')
plot3(vis_t2(1), vis_t2(2), vis_t2(3), 'bx')
for i=1:N_lines
    % Lines in camera 1
    plot3([vis_line_2D_1(1,i) vis_line_2D_1(4,i)], ...
          [vis_line_2D_1(2,i) vis_line_2D_1(5,i)], ...
          [vis_line_2D_1(3,i) vis_line_2D_1(6,i)], 'r')
    % Lines in camera 2
    plot3([vis_line_2D_2(1,i) vis_line_2D_2(4,i)], ...
          [vis_line_2D_2(2,i) vis_line_2D_2(5,i)], ...
          [vis_line_2D_2(3,i) vis_line_2D_2(6,i)], 'b')
    % Reference 3D lines
    plot3([line_3D(1,i) line_3D(4,i)], ...
          [line_3D(2,i) line_3D(5,i)], ...
          [line_3D(3,i) line_3D(6,i)], 'g')
    % Visualize normals of camera 1  
%     plot3([vis_t1(1) vis_n1(1,i)], ...
%           [vis_t1(2) vis_n1(2,i)], ...
%           [vis_t1(3) vis_n1(3,i)], 'k', 'LineWidth', 4)
    % Visualize normals of camera 2
%     plot3([vis_t2(1) vis_n2(1,i)], ...
%           [vis_t2(2) vis_n2(2,i)], ...
%           [vis_t2(3) vis_n2(3,i)], 'm', 'LineWidth', 4)
    % Visualize rh of each 3D line
%     plot3([-rh(1,i) rh(1,i)], ...
%           [-rh(2,i) rh(2,i)], ...
%           [-rh(3,i) rh(3,i)], '-.r')
    % Visualize m of each 3D line
%     plot3([-m(1,i) m(1,i)], ...
%           [-m(2,i) m(2,i)], ...
%           [-m(3,i) m(3,i)], '-.b')
      
    plot3([P1(1,i) P2(1,i)], ...
          [P1(2,i) P2(2,i)], ...
          [P1(3,i) P2(3,i)], 'c', 'LineWidth', 4)
    % Reconstructed 3D lines based on 3D points
%     plot3([P1(1,i) P2(1,i)], ...
%           [P1(2,i) P2(2,i)], ...
%           [P1(3,i) P2(3,i)], 'c', 'LineWidth', 1)
      
    % Reconstructed 3D lines based on 2D lines (this is what we need to make work)
    plot3([Pr1(1,i) Pr2(1,i)], ...
          [Pr1(2,i) Pr2(2,i)], ...
          [Pr1(3,i) Pr2(3,i)], 'k', 'LineWidth', 1)      
%     plot3(P1(1,i), P1(2,i), P1(3,i), 'r*', 'Markersize', 5)
%                 
%     plot3(P2(1,i), P2(2,i), P2(3,i), 'b*', 'Markersize', 5)
%                
%     plot3(Pr1(1,i), Pr1(2,i), Pr1(3,i), 'r.', 'Markersize', 10)
%                 
%     plot3(Pr2(1,i), Pr2(2,i), Pr2(3,i), 'b.', 'Markersize', 10)
end
visualize_camera_pose(camera, vis_R1, vis_t1)
visualize_camera_pose(camera, vis_R2, vis_t2)

%% Rotation error (paper Metrics for 3D Rotations: Comparison and Analysis)
% rot_error_k = norm(eye(3) - (C2_R_C1*new_R2'),'fro');
rot_error = norm(eye(3) - (C2_R_C1*R2'),'fro');
% disp(['Kanatani rotation error: ' num2str(rot_error_k)])
disp(['Plain rotation error:    ' num2str(rot_error)])

%% Translation error
% trans_error_k = norm(C2_P_C1 - new_t2);
trans_error = acos( (C2_P_C1' * t2) / (norm(C2_P_C1) * norm(t2)) );
% disp(['Kanatani translation error: ' num2str(trans_error_k)])
disp(['Plain translation error:    ' num2str(trans_error)])
