%% Generic script for 3D reconstruction from 2 views
% % Dimitris Zermas, 12/10/2015
% % PhD Candidate, Center for Distributed Robotics University of Minnesota
% % -----------------------------------------------------------------------
clearvars
clc

%% Select the two images --------------------------------------------------
[filename1, path_name1] = uigetfile('*.jpg','Select 1st Image');
[filename2, path_name2] = uigetfile('*.jpg','Select 2nd Image');

I1 = strcat(path_name1,filename1);
I2 = strcat(path_name2,filename2);

image1 = imread(I1);
image2 = imread(I2);

%% Select where to store the point cloud ----------------------------------
[out_path_name] = uigetdir('','Store reconstruction here');

%% Select used camera. This operation selects the apropriate intrinsics ---
% For the Olympus TG-4                       - 'Olympus'
% For reconstruction of trucks               - 'Trucks'
% For reconstruction of the Furukawa dataset - 'Furukawa'
camera = 'Trucks';

if strcmp(camera, 'Trucks')
    cam1_number = str2num(I1(end-5:end-4));
    cam2_number = str2num(I2(end-5:end-4));
else
    cam1_number = 0;
    cam2_number = 0;
end

%% Copmute the ASIFT corresponding points ---------------------------------
matchASIFT = ASIFT(image1, image2);

C1_pix = [matchASIFT(2:end,1)'; matchASIFT(2:end,2)'];
C2_pix = [matchASIFT(2:end,3)'; matchASIFT(2:end,4)'];
% Homogeneous pixels
C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];

%% Visualization of matches------------------------------------------------
% Show all
% match_plot(image1, image2, C1_pix', C2_pix');
% Visualize ASIFT points
% figure, imshow(image1), hold on;
% plot(C1_pix(1,:),C1_pix(2,:),'ro'), hold off
% figure, imshow(image2), hold on;
% plot(C2_pix(1,:),C2_pix(2,:),'bo'), hold off

%% Intrinsics--------------------------------------------------------------
[camera1, camera2] = camera_parameters_selection(camera, cam1_number, cam2_number);

%% Normalize points--------------------------------------------------------
% C*_normalized == 3 homogeneous vector
C1_n = normalize(C1_pix,camera1.fc,camera1.cc,camera1.kc,camera1.alpha_c);
C2_n = normalize(C2_pix,camera2.fc,camera2.cc,camera2.kc,camera2.alpha_c);
% C*_homogeneous == 2D in 3Dimensions, not normalized aka norm != 1
C1_h = [C1_n; ones(1,size(C1_n,2))];
C2_h = [C2_n; ones(1,size(C2_n,2))];

%% Motion parameters and depth estimation
% r1: 3D points with Kanatani correction
% r2: 3D points without Kanatani correction
[r1, r2, motion_params] = reconstruction(C1_h, C2_h);

h1 = motion_params.h1;
h2 = motion_params.h2;

%% Write to ply------------------------------------------------------------
pts2write = ['element vertex ' num2str(size(r1,2)) ' # number of points'];
plyname = [I1(end-5:end-4) '_' I2(end-5:end-4)];
fileIDw = fopen([out_path_name '/' filename1(1:end-6) plyname '.ply'],'w');
fprintf(fileIDw, '%s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n',...
    'ply', ...
    'format ascii 1.0', ...
    pts2write, ...
    'property float x', ...
    'property float y', ...
    'property float z', ...
    'property uchar diffuse_red', ...
    'property uchar diffuse_green', ...
    'property uchar diffuse_blue', ...
    'end_header');
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h1(1),h1(2),h1(3), 255,0,0);
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h2(1),h2(2),h2(3), 0,255,0);

for i=1:size(r1,2)
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', r1(1,i),r1(2,i),r1(3,i), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),1), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),2), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),3));
end
fclose(fileIDw);

system('rm -rf tmpASIFT*');
system('rm -rf tr*.txt');
system('rm -rf imgOut*');
system('rm -rf keys*');
system('rm -rf matchings*');