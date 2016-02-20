function [camera1, camera2] = camera_parameters_selection(str, cam1_number, cam2_number)

switch str
%% Camera Olympus Stylus TG-4
case 'Olympus'
    % Camera 1
    camera1.fc = [ 1422.84271   1425.39961 ];
    camera1.cc = [ 950.33846   737.83162 ];
    camera1.kc = [ 0.02694   -0.03312   0.00563   -0.00268  0.00000 ]; 
    camera1.alpha_c = 0.00000;
    % Camera 2
    camera2.fc = [ 1422.84271   1425.39961 ];
    camera2.cc = [ 950.33846   737.83162 ];
    camera2.kc = [ 0.02694   -0.03312   0.00563   -0.00268  0.00000 ]; 
    camera2.alpha_c = 0.00000;

%% Camera for the Furukawa dataset   
case 'Furukawa'
    % Camera 1
    camera1.fc = [ 6704.926882; 6705.241311 ]; % Focal length
    camera1.kc = [-0.125368; -0.097388; -0.003711; -0.000161; 0.000000]; % Non-linear distortions r = 1 + kc(1)r^2 + kc(2)r^4 + kc(3)r^6
    camera1.cc = [ 1038.251932; 957.560286 ]; % Principal point coordinates (cx,cy)
    camera1.alpha_c = 0.000103; % Skewness

    % Camera 2
    camera2.fc = [ 6682.125964; 6681.475962 ];
    camera2.kc = [-0.106090; -0.533543; -0.005174; 0.000517; 0.000000];
    camera2.cc = [ 1175.207200; 857.700292 ];
    camera2.alpha_c = 0.000101;

%% Cameras for the Truck Parking project
case 'Trucks'
    % Camera 1
    cams = readBundler('bundle.out');
    camera1.fc = [cams(cam1_number).focal cams(cam1_number).focal];
    camera1.kc = [cams(cam1_number).radial(1) cams(cam1_number).radial(2) 0 0 0];
    camera1.cc = [ 960; 540 ];
    camera1.alpha_c = 0;
    
    % Camera 2
    camera2.fc = [cams(cam2_number).focal cams(cam2_number).focal];
    camera2.kc = [cams(cam2_number).radial(1) cams(cam2_number).radial(2) 0 0 0];
    camera2.cc = [ 960; 540 ];
    camera2.alpha_c = 0;
end