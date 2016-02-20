% Import images
image1 = imread('../../data/trucks/truck01.jpg');
image2 = imread('../../data/trucks/truck02.jpg');

% image1 = imread('../../data/Furukawa_datasets/visualize/06.jpg');
% image2 = imread('../../data/Furukawa_datasets/visualize/14.jpg');

% Grayscale as needed by SIFT
I1 = single(rgb2gray(image1));
I2 = single(rgb2gray(image2));

%% Call the C code for the LSD!
% minimum length of line in pixels!
threshold = 10;

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

% ATTENTION! pts1 and pts2 have one-to-one correspondence!
C1_hom_pix = [F1(1,MATCHES(1,:)); F1(2,MATCHES(1,:)); ones(1, length(MATCHES))];
C2_hom_pix = [F2(1,MATCHES(2,:)); F2(2,MATCHES(2,:)); ones(1, length(MATCHES))];
clear MATCHES

% Visualize correspondences 1-by-1
for i=1:size(C1_hom_pix,2)
    figure, imshow([image1 image2]), hold on
    plot(C1_hom_pix(1,i), C1_hom_pix(2,i),'r*')
    plot(size(image1,2) + C2_hom_pix(1,i), C2_hom_pix(2,i), 'b*')
    pause, close all
end

% Visualize ASIFT points
figure, imshow(image1), hold on;
plot(C1_pix(1,:),C1_pix(2,:),'ro'), hold off
figure, imshow(image2), hold on;
plot(C2_pix(1,:),C2_pix(2,:),'bo'), hold off


%% Find corresponding lines
% Image 1
l = 1;
for i=1:size(lines1,1)
    % find the points closest to the i_th line in image1
    r = zeros(size(C1_hom_pix,2),1);
    for j=1:size(C1_hom_pix,2)
        r(j) = distPointToLineSegment( [lines1(i,1) lines1(i,2)], [lines1(i,3) lines1(i,4)], C1_hom_pix(1:2,j)' );
    end
    [dist, corr_ind] = min(r);
    % Keep only lines that are close to a corresponding point!
    if dist < 1
        corr_point1(l) = corr_ind;
        line_number1(l) = i;
        l = l+1;
        % Visualize line and closest point
%         visualise_LSD_lines(lines1(i,:), image1);
%         plot(C1_hom_pix(1,corr_ind),C2_hom_pix(2,corr_ind),'bo')
%         pause, close all
    end
end

% Image 2
l = 1;
for i=1:size(lines2,1)
    % find the points closest to the i_th line in image1
    r = zeros(size(C2_hom_pix,2),1);
    for j=1:size(C2_hom_pix,2)
        r(j) = distPointToLineSegment( [lines2(i,1) lines2(i,2)], [lines2(i,3) lines2(i,4)], C2_hom_pix(1:2,j)' );
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
        line_matches(count, 1:2) = [i id];
        figure, imshow([image1; image2]), hold on
        plot([final_lines1(i,1) final_lines1(i,3)], [final_lines1(i,2) final_lines1(i,4)], 'LineWidth',2,'Color','green')
        plot([final_lines2(id,1) final_lines2(id,3)], size(image1,1) + [final_lines2(id,2) final_lines2(id,4)], 'LineWidth',2,'Color','yellow')
        pause, close all
        count = count+1;
    end
end

%% These are rough estimates based on Ted's results
% Camera 1
% Focal length
fc = [1.9102798063e+03 1.9102798063e+03];
% Non-linear distortions r = 1 + kc(1)r^2 + kc(2)r^4 + kc(3)r^6
kc = [-1.3078984740e-01 1.0202590027e-01 0 0 0];
% Principal point coordinates (cx,cy) (assume center of camera)
cc = [size(image1,2)/2 size(image1,1)/2];
% Skewness (digital cameras usually are not skewed)
alpha_c = 1;

%% Generate more points on each line
npoints = 100;
figure(1), hold on;
figure(2), imshow(image1), hold on;
for i=1:size(line_matches,1)
    extra_pts1 = [extract_points_from_line(final_lines1(i,[1 3]), f_p1(i,:), npoints) ]';
    extra_pts2 = [extract_points_from_line(final_lines2(i,[1 3]), f_p2(i,:), npoints) ]';

    %% Normalise points--------------------------------------------------------
    npts1 = [normalize(extra_pts1,camera1.fc,camera1.cc,camera1.kc,camera1.alpha_c); ones(1,size(extra_pts1,2))];
    npts2 = [normalize(extra_pts2,camera2.fc,camera2.cc,camera2.kc,camera1.alpha_c); ones(1,size(extra_pts2,2))];
    
    %% Fundamental Matrix
%     [E, e1, e2] = fundmatrix(npts1, npts2);
% 
%     epipolar = diag(npts1'*E*npts2);
% 
%     %% Find Rotation and Translation {R,t} for the second camera---------------
%     [rot,t] = EssentialMatrixToCameraMatrix(E);
%     k=4;
%     vis_pts1 = []; vis_pts2 = [];
%     for k=1:4
%         R = rot(:,:,k);
%         h = t(:,:,k);
%         C1 = [0 0 0]';
%         C2 = h;
% 
%         for i=1:size(npts1,2)
%             vis_pts1(1:3,i) = npts1(:,i);
%             vis_pts2(1:3,i) = R*npts2(:,i) + h;
%         end
% 
%         plot3(C1(1),C1(2),C1(3),'rx')
%         xlabel('x'), ylabel('y'), zlabel('z')
%         axis equal, grid on, hold on
%         plot3(C2(1),C2(2),C2(3),'bx')
%         for i=1:5%size(v_k,2)
%             plot3([C1(1) vis_pts1(1,i)],[C1(2) vis_pts1(2,i)],[C1(3) vis_pts1(3,i)],'b')
%             plot3([C2(1) vis_pts2(1,i)],[C2(2) vis_pts2(2,i)],[C2(3) vis_pts2(3,i)],'g')
%         end
%         pause
%         close all
%     end
    
    %% Kanatani line reconstruction
    init_n1 = [f_p1(i,:) 1]';
    init_n2 = [f_p2(i,:) 1]';
    % Normalized to norm(n,2) = 1
    n1 = init_n1/norm(init_n1,2);

    n2 = init_n2/norm(init_n2,2);
    
    init_L = cross(n1, R*n2) / dot(h, R*n2);
    L(1:3,i) = init_L / norm(init_L,2);
    figure(1),plot3([L(1,i) 2*L(1,i)],[L(2,i) 2*L(2,i)],[L(3,i) 2*L(3,i)])
    figure(2),plot([final_lines1(i,1) final_lines1(i,3)],[final_lines1(i,2) final_lines1(i,4)],'LineWidth',2,'Color','yellow')
    pause
end

