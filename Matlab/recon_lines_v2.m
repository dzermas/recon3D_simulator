clearvars; clc;

% Import images
image1 = imread('../../data/trucks/truck01.jpg');
image2 = imread('../../data/trucks/truck02.jpg');

[out_path_name] = uigetdir('','Store reconstruction here');

% image1 = imread('../../Furukawa_datasets/visualize/06.jpg');
% image2 = imread('../../Furukawa_datasets/visualize/14.jpg');

% Grayscale as needed by SIFT
I1 = single(rgb2gray(image1));
I2 = single(rgb2gray(image2));

%% Call the C code for the LSD!
% minimum length of line in pixels!
threshold = 50;

[lines1,p1] = call_LSD('truck01', threshold);
% visualise_LSD_lines(lines1, image1);

[lines2,p2] = call_LSD('truck02', threshold);
% visualise_LSD_lines(lines2, image2);

%% Compute SIFT features near the line segment ends!
fc1 = [lines1(:,1)' lines1(:,3)'; lines1(:,2)' lines1(:,4)'; ones(1,2*length(lines1))*10;zeros(1,2*length(lines1))];
[F1,D1] = vl_sift(I1,'frames',fc1,'orientations');

fc2 = [lines2(:,1)' lines2(:,3)'; lines2(:,2)' lines2(:,4)'; ones(1,2*length(lines2))*10;zeros(1,2*length(lines2))];
[F2,D2] = vl_sift(I2,'frames',fc2,'orientations');

% Find correspondences
[MATCHES,SCORES] = vl_ubcmatch(D1, D2, 1.5);

% ATTENTION! C1_pix and C2_pix already have one-to-one correspondence!
C1_pix = [F1(1,MATCHES(1,:)); F1(2,MATCHES(1,:))];
C2_pix = [F2(1,MATCHES(2,:)); F2(2,MATCHES(2,:))];

% Visualize correspondences 1-by-1
% for i=1:size(C1_pix,2)
%     figure, imshow([image1 image2]), hold on
%     plot(C1_pix(1,i), C1_pix(2,i),'r*')
%     plot(size(image1,2) + C2_pix(1,i), C2_pix(2,i), 'b*')
%     pause, close all
% end

%% Intrinsics--------------------------------------------------------------
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

%% Estimate Motion Parameters
R1 = eye(3);
t1 = [0;0;0];
% Use the OpenGV library 5 point Niester algorithm
R_t = opengv ( 'fivept_nister_ransac', C1_b, C2_b );
R2 = R_t(:,1:3);
t2 = R_t(:,4);
% Compute Essential matrix
E = skew_sym(t2) * R2;

%% Find corresponding lines
% Image 1
l = 1;
for i=1:size(lines1,1)
    % find the points closest to the i_th line in image1
    r = zeros(size(C1_pix,2),1);
    for j=1:size(C1_pix,2)
        r(j) = distPointToLineSegment( [lines1(i,1) lines1(i,2)], [lines1(i,3) lines1(i,4)], C1_pix(:,j)' );
    end
    [dist, corr_ind] = min(r);
    % Keep only lines that are close to a corresponding point!
    if dist < 1
        corr_point1(l) = corr_ind;
        line_number1(l) = i;
        l = l+1;
        % Visualize line and closest point
        %         visualise_LSD_lines(lines1(i,:), image1);
        %         plot(pts1(1,corr_ind),pts1(2,corr_ind),'bo')
        %         pause, close all
    end
end

% Image 2
l = 1;
for i=1:size(lines2,1)
    % find the points closest to the i_th line in image2
    r = zeros(size(C2_pix,2),1);
    for j=1:size(C2_pix,2)
        r(j) = distPointToLineSegment( [lines2(i,1) lines2(i,2)], [lines2(i,3) lines2(i,4)], C2_pix(:,j)' );
    end
    [dist, corr_ind] = min(r);
    % Keep only lines that are close to a corresponding point!
    if dist < 1
        corr_point2(l) = corr_ind;
        line_number2(l) = i;
        l = l+1;
        % Visualize line and closest point
        %         visualise_LSD_lines(lines2(i,:), image2);
        %         plot(pts2(1,corr_ind),pts2(2,corr_ind),'bo')
        %         pause, close all
    end
end

% So we finalize which lines we are going to consider for the matching
final_lines1 = lines1(line_number1,:);
final_lines2 = lines2(line_number2,:);
f_p1 = p1(line_number1,:);
f_p2 = p2(line_number2,:);

% Match lines and visualize them
count = 1;
for i=1:size(corr_point1,2)
    [matches, id] = find(corr_point1(i) == corr_point2);
    if ~isempty(matches)
        figure, imshow([image1; image2]), hold on
        plot([final_lines1(i,1) final_lines1(i,3)], [final_lines1(i,2) final_lines1(i,4)], 'LineWidth',2,'Color','green')
        plot([final_lines2(id,1) final_lines2(id,3)], size(image1,1) + [final_lines2(id,2) final_lines2(id,4)], 'LineWidth',2,'Color','yellow')
        s = input('Keep line pair?', 's');
        if strcmp(s,'1')
            line_matches(count, 1:2) = [i id];
            count = count+1;
        end
        pause
        close all
    end
end

N_lines = size(line_matches, 1);

for i=1:N_lines
    L11_pix = [final_lines1(line_matches(i,1), 1) final_lines1(line_matches(i,1), 2)];
    L12_pix = [final_lines1(line_matches(i,1), 3) final_lines1(line_matches(i,1), 4)];
    L21_pix = [final_lines2(line_matches(i,2), 1) final_lines2(line_matches(i,2), 2)];
    L22_pix = [final_lines2(line_matches(i,2), 3) final_lines2(line_matches(i,2), 4)];
    
    L11_h = [(L11_pix - camera1.cc) ./ camera1.fc 1]';
    L12_h = [(L12_pix - camera1.cc) ./ camera1.fc 1]';
    L21_h = [(L21_pix - camera2.cc) ./ camera2.fc 1]';
    L22_h = [(L22_pix - camera2.cc) ./ camera2.fc 1]';
    % Pair up the lines between the two image plains
    line_2D_1(1:6,i) = [L11_h; L12_h];
    line_2D_2(1:6,i) = [L21_h; L22_h];
end

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
    Z1(i) = dot(cross(t2,R2*line_2D_2(1:3,i)), cross(line_2D_1(1:3,i),R2*line_2D_2(1:3,i))) / ...
        norm(cross(line_2D_1(1:3,i),R2*line_2D_2(1:3,i)))^2;
    
    Z2(i) = dot(cross(t2,R2*line_2D_2(4:6,i)), cross(line_2D_1(4:6,i),R2*line_2D_2(4:6,i))) / ...
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
    plot3([vis_line_2D_1(1,i) vis_line_2D_1(4,i)], ...
        [vis_line_2D_1(2,i) vis_line_2D_1(5,i)], ...
        [vis_line_2D_1(3,i) vis_line_2D_1(6,i)], 'r')
    
    plot3([vis_line_2D_2(1,i) vis_line_2D_2(4,i)], ...
        [vis_line_2D_2(2,i) vis_line_2D_2(5,i)], ...
        [vis_line_2D_2(3,i) vis_line_2D_2(6,i)], 'b')
    
    %     plot3([line_3D(1,i) line_3D(4,i)], ...
    %           [line_3D(2,i) line_3D(5,i)], ...
    %           [line_3D(3,i) line_3D(6,i)], 'g')
    %
    %     plot3([vis_t1(1) vis_n1(1,i)], ...
    %           [vis_t1(2) vis_n1(2,i)], ...
    %           [vis_t1(3) vis_n1(3,i)], 'k', 'LineWidth', 4)
    %
    %     plot3([vis_t2(1) vis_n2(1,i)], ...
    %           [vis_t2(2) vis_n2(2,i)], ...
    %           [vis_t2(3) vis_n2(3,i)], 'm', 'LineWidth', 4)
    %
    %     plot3([-rh(1,i) rh(1,i)], ...
    %           [-rh(2,i) rh(2,i)], ...
    %           [-rh(3,i) rh(3,i)], '-.r')
    %
    %     plot3([-m(1,i) m(1,i)], ...
    %           [-m(2,i) m(2,i)], ...
    %           [-m(3,i) m(3,i)], '-.b')
    
    %                 plot3([L1(1,i) L2(1,i)], ...
    %                     [L1(2,i) L2(2,i)], ...
    %                     [L1(3,i) L2(3,i)], 'c', 'LineWidth', 4)
    
    plot3([P1(1,i) P2(1,i)], ...
        [P1(2,i) P2(2,i)], ...
        [P1(3,i) P2(3,i)], 'c', 'LineWidth', 4)
    
    plot3([Pr1(1,i) Pr2(1,i)], ...
        [Pr1(2,i) Pr2(2,i)], ...
        [Pr1(3,i) Pr2(3,i)], 'k', 'LineWidth', 4)
    %
    %                 plot3(P1(1,i), P1(2,i), P1(3,i), 'b*', 'Markersize', 5)
    %
    %                 plot3(P2(1,i), P2(2,i), P2(3,i), 'b*', 'Markersize', 5)
    %
    %                 plot3(Pr1(1,i), Pr1(2,i), Pr1(3,i), 'r.', 'Markersize', 10)
    %
    %                 plot3(Pr2(1,i), Pr2(2,i), Pr2(3,i), 'r.', 'Markersize', 10)
end
visualize_camera_pose(camera1, vis_R1, vis_t1)
visualize_camera_pose(camera2, vis_R2, vis_t2)



pts2write = ['element vertex ' num2str(2*N_lines) ' # number of points'];
plyname = '01_02';
fileIDw = fopen([out_path_name '/lines_' plyname '.ply'],'w');
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
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', t1(1),t1(2),t1(3), 255,0,0);
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', t2(1),t2(2),t2(3), 0,255,0);

for i=1:N_lines
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', Pr1(1,i),Pr1(2,i),Pr1(3,i), 255,255,0);
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', Pr2(1,i),Pr2(2,i),Pr2(3,i), 255,255,0);
%     fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', P1(1,i),P1(2,i),P1(3,i), 255,0,255);
%     fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', P2(1,i),P2(2,i),P2(3,i), 255,0,255);
end
fclose(fileIDw);

system('rm -rf tmpASIFT*');
system('rm -rf imgOut*');
system('rm -rf keys*');
system('rm -rf matchings*');
