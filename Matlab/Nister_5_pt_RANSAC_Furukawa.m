%% Nister 5 point with RANSAC outputs the motion parameters in a P = [R|t] format
clearvars
clc

%% Select corresponding points between two images -------------------------
I1 = '../../data/Furukawa_datasets/visualize/06.jpg';
I2 = '../../data/Furukawa_datasets/visualize/14.jpg';

[out_path_name] = uigetdir('','Store reconstruction here');

image1 = imread(I1);
image2 = imread(I2);

%% Import the ASIFT results -----------------------------------------------
matchASIFT = ASIFT(image1, image2);

C1_pix = [matchASIFT(2:end,1)'; matchASIFT(2:end,2)'];
C2_pix = [matchASIFT(2:end,3)'; matchASIFT(2:end,4)'];

C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];


%% Intrinsics--------------------------------------------------------------
% % Camera 1
% camera1.fc = [ 6704.926882; 6705.241311 ]; % Focal length
% camera1.kc = [-0.125368; -0.097388; -0.003711; -0.000161; 0.000000]; % Non-linear distortions r = 1 + kc(1)r^2 + kc(2)r^4 + kc(3)r^6
% camera1.cc = [ 1038.251932; 957.560286 ]; % Principal point coordinates (cx,cy)
% camera1.alpha_c = 0.000103; % Skewness
% 
% % Camera 2
% camera2.fc = [ 6682.125964; 6681.475962 ];
% camera2.kc = [-0.106090; -0.533543; -0.005174; 0.000517; 0.000000];
% camera2.cc = [ 1175.207200; 857.700292 ];
% camera2.alpha_c = 0.000101;

%% Play around with intrinsics to see how they affect result
% Camera 1
camera1.fc = [ 6705; 6705 ]; % Focal length
camera1.kc = [-0.125368; -0.097388; 0; 0; 0.000000]; % Non-linear distortions r = 1 + kc(1)r^2 + kc(2)r^4 + kc(3)r^6
camera1.cc = [ 975; 875 ]; % Principal point coordinates (cx,cy)
camera1.alpha_c = 0; % Skewness

% Camera 2
camera2.fc = [ 6682; 6682 ];
camera2.kc = [-0.106090; -0.533543; 0; 0; 0.000000];
camera2.cc = [ 975; 875 ];
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
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );

%% Find Rotation and Translation {R,t} for the second camera---------------
R2 = R_t(:,1:3);
h2 = R_t(:,4);
Essential = skew_sym(h2) * R2; % E = skew_sym(t) * R
% [rot,t] = EssentialMatrixToCameraMatrix(Essential);

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

%% Consider solution only if it is reconstructing points in front of the two cameras
% for i=1:size(v_k,2)
%     vis_pts1(1:3,i) = R1*v_k(:,i) + h1;
%     vis_pts2(1:3,i) = R2*u_k(:,i) + h2;
% end
% 
% figure, xlabel('x'), ylabel('y'), zlabel('z')
% axis equal, grid on, hold on
% plot3(h1(1),h1(2),h1(3),'rx')
% plot3(h2(1),h2(2),h2(3),'bx')
% for i=1:10:size(v_k,2)
%     plot3([h1(1) vis_pts1(1,i)],[h1(2) vis_pts1(2,i)],[h1(3) vis_pts1(3,i)],'b')
%     plot3([h2(1) vis_pts2(1,i)],[h2(2) vis_pts2(2,i)],[h2(3) vis_pts2(3,i)],'g')
%     plot3(r1(1,i),r1(2,i),r1(3,i),'r.')
%     plot3(r2(1,i),r2(2,i),r2(3,i),'c.')
%     plot3(r3(1,i),r3(2,i),r3(3,i),'b.')
% end
% pause

%% Write to ply------------------------------------------------------------
pts2write = ['element vertex ' num2str(size(r1_filtered,2)) ' # number of points'];
plyname = [I1(end-5:end-4) '_' I2(end-5:end-4)];
fileIDw = fopen([out_path_name '/predator_' plyname '.ply'],'w');
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