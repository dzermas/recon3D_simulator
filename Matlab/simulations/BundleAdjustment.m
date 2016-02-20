%% Bundle Adjustment
clearvars
clc

%% Initialize
depth_range = 10;
N_points = 5;
N_cameras = 10;
focal_length = 600;

% Generate 3D points
for np=1:N_points
    P{1}(1:3,np) = [(2 * rand)-1 (2 * rand)-1 ((2 * rand)-1)]';
    R{1} = eye(3);
    t{1} = zeros(3,1);
end

for nc=2:N_cameras
    % Generate Rotations
    R{nc} = generate_random_rotation;
    % Generate Translations
    t{nc} = [(2 * rand)-1; (2 * rand)-1 ; (2 * rand)-1];
end

% Generate 3D points in other cameras
for nc=1:N_cameras
    for np=1:N_points
        noise = 1/focal_length * [rand rand 0]';
        P{nc}(1:3,np) = R{nc}*P{1}(:,np) + t{nc};
        u{nc}(1:3,np) = P{nc}(:,np) ./ P{nc}(3,np) + noise;
    end
end

% Dummy initialization of BA
% Find depth estimation for the 3D points wrt camera 1
for np=1:N_points
    A = [];
    b = [];
    for nc=2:N_cameras
        A = [A ; skew_sym(u{nc}(:,np)) * R{nc} * u{1}(:,np)];
        b = [b ; - skew_sym(u{nc}(:,np)) * t{nc}];
    end
    depth_est(np) = A\b;
end

%% BA implementation
% Initial solution
for np=1:N_points
    P0(1:3,np) = u{1}(:,np) * depth_est(np);
end

% Noisy homogeneous coordinate

while 

