%% Nister 5 point with RANSAC outputs the motion parameters in a P = [R|t] format
clear all
clc

%% Select corresponding points between two images -------------------------
I1 = '../../data/trucks/truck02.jpg';
I2 = '../../data/trucks/truck01.jpg';

[out_path_name] = uigetdir('','Store reconstruction here');

cam1_number = str2num(I1(end-5:end-4));
cam2_number = str2num(I2(end-5:end-4));

image1 = imread(I1);
image2 = imread(I2);

tic
%% Import the ASIFT results -----------------------------------------------
matchASIFT = ASIFT(image1, image2);

C1_pix = [matchASIFT(2:end,1)'; matchASIFT(2:end,2)'];
C2_pix = [matchASIFT(2:end,3)'; matchASIFT(2:end,4)'];

C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];

%% Truck cams intrinsics---------------------------------------------------
cams = readBundler('bundle.out');
camera1.fc = [cams(cam1_number).focal cams(cam1_number).focal];
camera1.kc = [cams(cam1_number).radial(1) cams(cam1_number).radial(2) 0 0 0];
camera1.cc = [ 960; 540 ];
camera1.alpha_c = 0;

camera2.fc = [cams(cam1_number).focal cams(cam1_number).focal];
camera2.kc = [cams(cam1_number).radial(1) cams(cam1_number).radial(2) 0 0 0];
camera2.cc = [ 960; 540 ];
camera2.alpha_c = 0;

%% Normalize points--------------------------------------------------------
% C*_normalized == 3 homogeneous vector
C1_n = normalize(C1_pix,camera1.fc,camera1.cc,camera1.kc,camera1.alpha_c);
C2_n = normalize(C2_pix,camera2.fc,camera2.cc,camera2.kc,camera2.alpha_c);
% C*_homogeneous == 2D in 3Dimensions, not normalized aka norm != 1
C1_h = [C1_n; ones(1,size(C1_n,2))];
C2_h = [C2_n; ones(1,size(C2_n,2))];
% C*_bearing == normalized 3 homogeneous vector
for i=1:size(C1_h,2)
    C1_b(1:3,i) = C1_h(:,i) / norm(C1_h(:,i));
    C2_b(1:3,i) = C2_h(:,i) / norm(C2_h(:,i));
end

%% Try OpenGV Kneip
tic
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );
toc
%% Find Rotation and Translation {R,t} for the second camera---------------
R2 = R_t(:,1:3);
h2 = R_t(:,4);

Essential = skew_sym(h2) * R2; % E = skew_sym(t) * R

%% Correct points to accurately validate the epipolar constraint-----------
[v_k, u_k, Dv_k, Du_k, keepers] = correct_points_kanatani(C1_h, C2_h, Essential);

%% Triangulate kanatani
R1 = eye(3);
h1 = [0 0 0]';

for i=1:size(v_k,2)
    % No Kanatani correction reconstruction
    Z1(i) = dot(cross(h2,R2*C2_h(:,i)), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
    Z2(i) = dot(cross(C1_h(:,i),h2), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
    
    % Kanatani correction reconstruction
    Zk1(i) = dot(cross(h2,R2*u_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
    
    Zk2(i) = dot(cross(h2,v_k(:,i)), cross(v_k(:,i),R2*u_k(:,i))) / ...
        norm(cross(v_k(:,i),R2*u_k(:,i)))^2;
           
    n = cross(h2,v_k(:,i))/norm(cross(h2,v_k(:,i)));
    m = cross(n,R2*u_k(:,i));
    
    DZk(i) = -(dot(m, Zk1(i)*Dv_k(:,i) - Zk2(i)*R2*Du_k(:,i))) / dot(m,v_k(:,i));
    
    r1(1:3,i) = (Zk1(i) + DZk(i))*v_k(:,i);
    r2(1:3,i) = Z1(i)*C1_h(:,i);
end

%% Reprojection error
error_k = reprojection_error(R1, h1, R2, h2, r1, C1_b(:,keepers), C2_b(:,keepers));
error   = reprojection_error(R1, h1, R2, h2, r2, C1_b(:,keepers), C2_b(:,keepers));
disp(['Kanatani reprojection error: ' num2str(error_k)])
disp(['Plain reprojection error:    ' num2str(error)])

%% Filter reconstructed points based on depth
for i=1:length(Z1)
    if Z1(i) < 2*mean(Z1)
        r1_filtered(1:3,i) = r1(:,i);
    end
end

toc

%% Write to ply------------------------------------------------------------
pts2write = ['element vertex ' num2str(size(r1_filtered,2)) ' # number of points'];
plyname = [I1(end-5:end-4) '_' I2(end-5:end-4)];
fileIDw = fopen([out_path_name '/trucks_' plyname '.ply'],'w');
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

for i=1:size(r1_filtered,2)
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', r1_filtered(1,i),r1_filtered(2,i),r1_filtered(3,i), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),1), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),2), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),3));
end
fclose(fileIDw);

system('rm -rf tmpASIFT*');
system('rm -rf imgOut*');
system('rm -rf keys*');
system('rm -rf matchings*');