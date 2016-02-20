clearvars; clc;

% Import images
image1 = imread('../../data/trucks/tr01.jpg');
image2 = imread('../../data/trucks/tr02.jpg');

% image1 = imread('../../Furukawa_datasets/visualize/06.jpg');
% image2 = imread('../../Furukawa_datasets/visualize/14.jpg');

% Grayscale as needed by SIFT
I1 = single(rgb2gray(image1));
I2 = single(rgb2gray(image2));

%% Call the C code for the LSD!
% minimum length of line in pixels!
threshold = 50;

[lines1,p1] = call_LSD('tr01', threshold);
% visualise_LSD_lines(lines1, image1);

[lines2,p2] = call_LSD('tr02', threshold);
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
% for i=1:size(C1_hom_pix,2)
%     figure, imshow([image1 image2]), hold on
%     plot(C1_pix(1,i), C1_pix(2,i),'r*')
%     plot(size(image1,2) + C2_pix(1,i), C2_pix(2,i), 'b*')
%     pause, close all
% end

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
figure, imshow([image1; image2]), hold on
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
        close all
    end
end

%% Intrinsics--------------------------------------------------------------
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
    % normals expressed in the world coordinates
    n1(:,i) = init_n1/norm(init_n1,2);
    n2(:,i) = init_n2/norm(init_n2,2);

    % {m,rh}-representation of the 3D line
    p(:,i) = cross(n1(:,i), R2*n2(:,i)) / dot(t2, R2*n2(:,i));
    m(:,i) = p(:,i) / norm(p(:,i));
    rh(:,i) = cross(p(:,i), n1(:,i)) / norm(p(:,i))^2;
    
    % Begin (L1) and end (L2) of the reconstructed line
    L1(:,i) = rh(:,i) - m(:,i);
    L2(:,i) = rh(:,i) + m(:,i);
end

