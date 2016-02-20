clearvars, clc

%% Point from 0 to 10
depth_range = 10;

camera.fc = [400;400];
camera.cc = [320;240];
camera.dx = 640;
camera.dy = 480;

N_points = 10;
C1_P = []; %% Euclidean
C1_pix = []; %% Pixel
C1_h = []; %% Homogeneous
C1_b = []; %% Bearings
C1_d = []; %% Depths

C2_P = []; %% Euclidean
C2_pix = []; %% Pixel
C2_h = []; %% Homogeneous
C2_b = []; %% Bearings
C2_d = []; %% Depths

%% Rotations
random_angle_1 = rand * (-5) * pi / 180;
random_angle_2 = rand * (-5) * pi / 180;
random_angle_3 = rand * (-5) * pi / 180;

R_Z = [cos(random_angle_3) -sin(random_angle_3) 0;
    sin(random_angle_3) cos(random_angle_3) 0;
    0 0 1];

R_Y = [cos(random_angle_2) 0 sin(random_angle_2);
        0 1 0;
        -sin(random_angle_2) 0 cos(random_angle_2)];

R_X = [1 0 0;
    0 cos(random_angle_1) -sin(random_angle_1);
    0 sin(random_angle_1) cos(random_angle_1)];

C2_R_C1 = R_Z;

%% Translations
C2_P_C1 = [0.5 * rand; 0 ; 0];

for i = 1:N_points
    random_x_pixel = round(rand * camera.dx);
    random_y_pixel = round(rand * camera.dy);
    random_depth = 1.5 + round(rand * depth_range);
    pixel = [random_x_pixel;random_y_pixel];
    homogeneous = [(pixel - camera.cc) ./ camera.fc ; 1];
    bearing = homogeneous / norm(homogeneous);
    euclidean = bearing * random_depth;
    C1_P = [C1_P euclidean];
    C1_d = [C1_d random_depth];
    C1_h = [C1_h homogeneous];
    C1_b = [C1_b bearing];
    C1_pix = [C1_pix pixel];
    
    C2_euclidean = C2_R_C1 * euclidean + C2_P_C1;
    C2_homogeneous = C2_euclidean / C2_euclidean(3);
    C2_bearing = C2_homogeneous / norm(C2_homogeneous);
    C2_pixel = C2_homogeneous(1:2) .* camera.fc + camera.cc;
    
    C2_P = [C2_P C2_euclidean];
    C2_h = [C2_h C2_homogeneous];
    C2_b = [C2_b C2_bearing];
    C2_pix = [C2_pix C2_pixel];
end

%% Essential_For_Bearing
E_bearings =  skew_sym(C2_P_C1) * C2_R_C1;
[Essential, ~, ~] = fundmatrix(C1_h, C2_h);

for i = 1:N_points
    assert(norm( C1_h(:,i)' * E_bearings' * C2_h(:,i)) < 1e-10)
end

%% Essential For homogeneous
for i = 1:N_points
    assert(norm( C1_h(:,i)' * Essential' * C2_h(:,i)) < 1e-10)
end

%% Fundamental
K = [camera.fc(1) 0 camera.cc(1);
    0 camera.fc(2) camera.cc(2);
    0 0 1];

F_matrix = (inv(K)') * Essential * inv(K);

%% Essential For Bearing
for i = 1:N_points
    assert(norm( [C1_pix(:,i);1]' * F_matrix' * [C2_pix(:,i);1]) < 1e-10)
end

% [rot,t] = EssentialMatrixToCameraMatrix(Essential);
% R2 = rot(:,:,1);
% h2 = t(:,:,1);

%% Kanatani {R,h}
[U,~,~] = svd(Essential*Essential');
h2 = U(:,3)./max(U(:,3));

K = -skew_sym(h2)*Essential;

[U,~,V] = svd(K);

R2 = U*diag([1 1 det(U*V')])*V';

R1 = eye(3);
h1 = [0 0 0]';
% R2 = C2_R_C1;
% h2 = C2_P_C1 / norm(C2_P_C1);


for i=1:size(C1_b,2)
    Z1(i) = -dot(cross(h2,R2*C2_b(:,i)), cross(C1_b(:,i),R2*C2_b(:,i))) / ...
        norm(cross(C1_b(:,i),R2*C2_b(:,i)))^2;
    
    Z2(i) = -dot(cross(h2,C1_b(:,i)), cross(C1_b(:,i),R2*C2_b(:,i))) / ...
        norm(cross(C1_b(:,i),R2*C2_b(:,i)))^2;
    
    r1(1:3,i) = Z1(i)*C1_b(:,i);
    r2(1:3,i) = Z2(i)*C2_b(:,i);
end

%% Express all points w.r.t. the first camera for visualization purposes
vis_pts1 = []; vis_pts2 = [];
for i=1:size(C1_b,2)
    vis_pts1(1:3,i) = C1_b(:,i);
    vis_pts2(1:3,i) = R2*C2_b(:,i) + h2;
    vis_r1(1:3,i) = r1(1:3,i);
    vis_r2(1:3,i) = r2(1:3,i);
end

figure, xlabel('x'), ylabel('y'), zlabel('z')
axis equal, grid on, hold on
plot3(h1(1), h1(2), h1(3), 'rx')
plot3(h2(1), h2(2), h2(3), 'bx')
for i=1:size(C1_b,2)
    plot3([h1(1) vis_pts1(1,i)], [h1(2) vis_pts1(2,i)], [h1(3) vis_pts1(3,i)], 'b')
    plot3([h2(1) vis_pts2(1,i)], [h2(2) vis_pts2(2,i)], [h2(3) vis_pts2(3,i)], 'g')
    plot3(C1_P(1,i), C1_P(2,i), C1_P(3,i), 'b^')
    plot3(vis_r1(1,i), vis_r1(2,i), vis_r1(3,i), 'b.')
    plot3(vis_r2(1,i), vis_r2(2,i), vis_r2(3,i), 'r.')
end