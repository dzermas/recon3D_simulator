%% #---- Educational platform for 3D reconstruction from two images ----#
% Dimitris Zermas, 2/15/2016
% University of Minnesota, Center for Distributed Robotics

% Feel free to explore and play with this code. This is for educational
% purposes only, hope you find it helpful!

clearvars
clc
close all

addpath('common')
%% Initialize
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

%% Camera Intrinsics
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
C2_R_C1 = generateRandomRotation();

% Translations
% reference t1
C1_P_C1 = [0; 0; 0];
% reference t2
C2_P_C1 = [2*rand; rand ; 0.2*rand];

F = skewSym(C2_P_C1) * C2_R_C1;

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

%% Visualize objects and cameras
f1 = figure;
hold on
plot3(xyz(1,:), xyz(2,:), xyz(3,:),'.')
visualizeCameraPose(camera, C2_R_C1, C2_P_C1, 'r')
visualizeCameraPose(camera, C1_R_C1, C1_P_C1, 'b') 
axis equal
title('Real 3D points')
hold off
movegui(f1,'southwest')

%% Generate points on the image planes
noise_gain = 0.001; % Try to see how much is considered "too much noise" :)
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
f2 = figure;
hold on
subplot(1,2,1);
plot(C1_pix(1,:), C1_pix(2,:), 'r.')
title('What the left camera sees')
subplot(1,2,2);
plot(C2_pix(1,:), C2_pix(2,:), 'b.')
title('What the right camera sees')
hold off
movegui(f2,'south')

[E, residual] = computeEssentialHL(C1_h, C2_h);
disp(['Residual computing the Essential matrix :' num2str(residual)]);

[R, t] = motionParametersKanatani(E, C1_h, C2_h);

[recon_kanatani, recon] = triangulationKanatani(C1_h, C2_h, R, t);

repr_err = reprojectionError(C1_R_C1, C1_P_C1, R, t, recon, C1_h, C2_h);
disp(['Rotation error :' num2str(norm(C2_P_C1' * R - eye(3), 'fro'))]);
disp(['Reprojection error :' num2str(repr_err)]);

% Plot 3D
f3 = figure;
hold on
title('Reconstruction')
plot3(C1_P(1,:), C1_P(2,:), C1_P(3,:), 'b.')
plot3(recon_kanatani(1,:), recon_kanatani(2,:), recon_kanatani(3,:), 'r*')
plot3(recon(1,:), recon(2,:), recon(3,:), 'g*')
legend({'real points', 'reconstructed with Kanatani correction', 'reconstructed points'})
axis equal
visualizeCameraPose(camera, C1_R_C1, C1_P_C1, [1 0 0]), text(C1_P_C1(1), ...
    C1_P_C1(2), C1_P_C1(3), 'C1ref', 'FontSize', 10, 'Color', 'k');
visualizeCameraPose(camera, C2_R_C1, C2_P_C1, [0 1 0]), text(C2_P_C1(1), ...
    C2_P_C1(2), C2_P_C1(3), 'C2ref', 'FontSize', 10, 'Color', 'k');
visualizeCameraPose(camera, R, t, [0 0 1]), text(t(1), t(2), t(3), 'C2', ...
    'FontSize', 10, 'Color', 'k');
hold off
movegui(f3,'southeast')
