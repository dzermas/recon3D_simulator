clearvars
clc
close all

addpath('common')
%% Initialize
dense = 0.5; % How dense the reconstruction should be
nClasses = 2; % number of SOMs to run
cube_points = 50;
sphere_points = 100;
N_points = cube_points + sphere_points;

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
camera.kc = [0 0 0 0 0];
camera.alpha_c = 0;
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
C2_P_C1 = [2*rand; rand ; 0.2*rand];

F = skew_sym(C2_P_C1) * C2_R_C1;

%% Create the 3D points
% Cube
cube_translation = [-1 0 7]; % x, y, z

x = linspace(-1, 1, sqrt(N_points));
y = linspace(-1, 1, sqrt(N_points));
z = -1;

[X, Y, Z] = ndgrid(x,y,z);
temp = [X(:),Y(:),Z(:)];
temp = cat(1,temp(:,[1 2 3]),temp(:,[2 3 1]));
temp = unique(temp,'rows');
DR = rotation ( 'x', 20 );
for i=1:length(temp)
    xyz(1:3,i) = DR*temp(i,:)' + cube_translation';
end

% Sphere
sphere_translation = [3 0 7]; % x, y, z

[X_temp, Y_temp, Z_temp] = sphere(30);
idx = find(Z_temp < -0.50);
X = X_temp(idx) + sphere_translation(1)*ones(length(idx),1);
Y = Y_temp(idx) + sphere_translation(2)*ones(length(idx),1);
Z = Z_temp(idx) + sphere_translation(3)*ones(length(idx),1);
xyz = [xyz [X Y Z]'];

%% Visualize object and cameras
figure, 
plot3(xyz(1,:), xyz(2,:), xyz(3,:),'.'), hold on
visualize_camera_pose(camera, C2_R_C1, C2_P_C1, 'r'),
visualize_camera_pose(camera, C1_R_C1, C1_P_C1, 'b'), axis equal, 
title('Initial 3D shapes')
hold off

%% Generate points on the image planes
noise_gain = 0.0;
for i = 1:size(xyz,2)
    euclidean = xyz(:,i);
    homogeneous = euclidean / euclidean(3);
    bearing = homogeneous / norm(homogeneous);
    pixel = homogeneous(1:2).* camera.fc + camera.cc;
    
    % Noise in measurements 1
    g_noise_1 = [noise_gain*rand; noise_gain*rand; 0];
    
    homogeneous_noisy = homogeneous + g_noise_1;
    bearing_noisy = homogeneous_noisy / norm(homogeneous_noisy);
    
    C1_P(:,i) = euclidean;
    C1_h(:,i) = homogeneous_noisy;
    C1_b(:,i) = bearing_noisy;
    C1_pix(:,i) = pixel(1:2);
    
    C2_euclidean = C2_R_C1' * (euclidean) - C2_R_C1 * C2_P_C1;
    C2_homogeneous = C2_euclidean / C2_euclidean(3);
    C2_bearing = C2_homogeneous / norm(C2_homogeneous);
    C2_pixel = C2_homogeneous(1:2) .* camera.fc + camera.cc;
    
    C2_P(:,i) = C2_euclidean;
    C2_h(:,i) = C2_homogeneous;
    C2_b(:,i) = C2_bearing;
    C2_pix(:,i) = C2_pixel;
end

%% Visualize what the cameras see!!!
figure, hold on
subplot(1,2,1);
plot(C1_pix(1,:), C1_pix(2,:), 'ro')
title('What the left camera sees')
subplot(1,2,2);
plot(C2_pix(1,:), C2_pix(2,:), 'bo')
title('What the right camera sees')
hold off

E = ComputeEssentialHL(C1_h, C2_h);
[R, t, error] = motion_parameters(E, C1_h, C2_h);

recon = reconstructionSimple(C1_h, C2_h, R, t);

error = reprojection_error(C1_R_C1, C1_P_C1, R, t, recon, C1_h, C2_h);
disp(['Rotation error :' num2str(norm(C2_P_C1' * R - eye(3), 'fro'))]);

% Plot 3D
figure, hold on
title('Reconstruction')
plot3(C1_P(1,:), C1_P(2,:), C1_P(3,:), 'b.')
plot3(recon(1,:), recon(2,:), recon(3,:), 'r*'), axis equal
% plot3(r2(1,:), r2(2,:), r2(3,:), 'go')
visualize_camera_pose(camera, C1_R_C1, C1_P_C1, [1 0 0]), text(C1_P_C1(1), C1_P_C1(2), C1_P_C1(3), 'C1ref', 'FontSize', 10, 'Color', 'k');
visualize_camera_pose(camera, C2_R_C1, C2_P_C1, [0 1 0]), text(C2_P_C1(1), C2_P_C1(2), C2_P_C1(3), 'C2ref', 'FontSize', 10, 'Color', 'k');
visualize_camera_pose(camera, R, t, [0 0 1]), text(t(1), t(2), t(3), 'C2', 'FontSize', 10, 'Color', 'k');
hold off