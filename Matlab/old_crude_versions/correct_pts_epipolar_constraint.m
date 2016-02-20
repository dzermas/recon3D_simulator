%% Select corresponding points between two images -------------------------
image1 = imread('../../data/adel_left_1070x750.jpg');
image2 = imread('../../data/adel_right_1070x750.jpg');

%% Import the ASIFT results -----------------------------------------------
load('../../data/match_ASIFT.mat');
C1_pix = [matchASIFT(:,1)'; matchASIFT(:,2)'];
C2_pix = [matchASIFT(:,3)'; matchASIFT(:,4)'];

C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];

%% OR Manually select points-----------------------------------------------
% Alternatively select points from the two images starting from the first
% and ending on the second
% eg I1 -> I2 -> I1 -> I2 ...
% figure; imshow([image1; image2])
% [x1, y1] = ginput;
% 
% C1_pix = [x1(1:2:end)'; y1(1:2:end)'];
% C2_pix = [x1(2:2:end)'; y1(2:2:end)' - size(image1,1)];
% 
% C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
% C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];

% Already picked points ;)
% load('../../data/manually_selected_pts2.mat')


%% OR Select points through SIFT-------------------------------------------
% I1 = single(rgb2gray(image1));
% I2 = single(rgb2gray(image2));
% % Compute SIFT features
% [F1,D1] = vl_sift(I1);
% [F2,D2] = vl_sift(I2);
% % Find correspondences
% [MATCHES,SCORES] = vl_ubcmatch(D1, D2, 4.5);
% % Take all the correspondences to see the results
% corrs = MATCHES;
% C1_pix = F1(1:2,corrs(1,:));
% C2_pix = F2(1:2,corrs(2,:));
% 
% C1_hom_pix = [C1_pix; ones(1,length(C1_pix))];
% C2_hom_pix = [C2_pix; ones(1,length(C2_pix))];

%% Visualization of matches------------------------------------------------
% Show all
% match_plot(image1, image2, C1_pix', C2_pix');

%% Intrinsics--------------------------------------------------------------
% These are rough estimates based on Ted's results
% Camera 1
camera1.fc = [1.9102798063e+03 1.9102798063e+03]; % Focal length
camera1.kc = [-1.3078984740e-01 7.7652638624e-02 0 0 0]; % Non-linear distortions r = 1 + kc(1)r^2 + kc(2)r^4 + kc(3)r^6
camera1.cc = [size(image1,2)/2 size(image1,1)/2]; % Principal point coordinates (cx,cy) (assume center of camera)
camera1.alpha_c = 0; % Skewness (digital cameras usually are not skewed)

% Camera 2
camera2.fc = [1.9180646963e+03 1.9180646963e+03];
camera2.kc = [-1.2658704927e-01 5.7129056080e-02 0 0 0];
camera2.cc = [size(image1,2)/2 size(image1,1)/2];
camera2.alpha_c = 0;

% Generic camera model
camera.fc = [10 10];
camera.kc = [0 0 0 0 0];
camera.cc = [size(image1,2)/2 size(image1,1)/2];
camera.alpha_c = 0;

%% Olympus intrinsics
load('../../olympus_calibrate/Olympus_intrinsics.mat');

%% Normalise points--------------------------------------------------------
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

%% Essential Matrix -------------------------------------------------------
% For normalized points, fundamental matrix is essential matrix
[Essential, e1, e2] = fundmatrix(C1_h, C2_h);
[Essential_pix, e1_pix, e2_pix] = fundmatrix(C1_hom_pix, C2_hom_pix);

% Is the epipolar constraint validated?
epipolar = diag(C2_h'*Essential*C1_h);
epipolar_pix = diag(C2_hom_pix'*Essential_pix*C1_hom_pix);

residual = epipolar - epipolar_pix;

%% Find Rotation and Translation {R,t} for the second camera---------------
[rot,t] = EssentialMatrixToCameraMatrix(Essential);
% Select the correct out of 4 possible solutions

%% Visualize the epipolar lines--------------------------------------------
% Need to express in pixel scaled homogeneous coordinates ++ Essential
% matrix should be computed from the homogeneous pixel points (C*_hom_pix)
% for i=1:size(C1_hom_pix,2)
%     line1 = C2_hom_pix(:,i)'*Essential_pix;
%     line2 = Essential_pix*C1_hom_pix(:,i);
%     % For image1
%     figure, imshow(image1), hold on
%     plot(C1_hom_pix(1,i),C1_hom_pix(2,i),'rx')
%     plot([0 2000], [-line1(3)/line1(2) -line1(1)*2000/line1(2)-line1(3)/line1(2)],'b'),
%     hold off
%     pause
%     % For image2
%     figure, imshow(image2), hold on
%     plot(C2_hom_pix(1,i),C2_hom_pix(2,i),'yx')
%     plot([0 2000], [-line2(3)/line2(2) -line2(1)*2000/line2(2)-line2(3)/line2(2)],'b')
%     hold off
%
%     pause
%
%     close all
% end


%% Correct points to accurately validate the epipolar constraint-----------
% Covariances show the variance of the points
% Assume independent noise in each direction on the image (x, y)
k = [0 0 1]';
Vv = eye(3) - k*k';
Vu = eye(3) - k*k';

% Iterate over all the points
% v_c: points in image 1
% u_c: points in image 2
J = []; v_c = []; u_c = [];
for i=1:length(C1_b)
    v = C1_h(:,i); u = C2_h(:,i);
    
    % Compute the error of the points in the two images
    Dv = (u'*Essential'*v * Vv * Essential * u) / (norm(Vv * Essential' * v, 2)^2 + norm(Vu * Essential * u,2)^2);
    Du = (u'*Essential'*v * Vu * Essential' * v) / (norm(Vv * Essential' * v, 2)^2 + norm(Vu * Essential * u,2)^2);
    
    % Compute the residual to see whether the points correspond or not
    v_c(1:3,i) = v - Dv; u_c(1:3,i) = u - Du;
    J(i) = (v'*Essential*u)^2 / (norm(Vv * Essential' * v_c(1:3,i), 2)^2 + norm(Vu * Essential * u_c(1:3,i),2)^2);
end

% Is the epipolar constraint validated with the corrected points?
epipolar_c = diag(v_c'*Essential*u_c);


%% Is kept with confidence 80% if J < CritVal (Degrees of Freedom = 1) ----
keep = J < 0.0642;
sum(keep)

v_k = v_c(:, keep);
u_k = u_c(:, keep);

% v_k = C1_h;
% u_k = C2_h;

%% Visualize the space lines and where they should meet--------------------
% Assume the first camera has center C1 = [0,0,0] and the second C2 = C1 + h
% where h is the translation of the 2nd camera w.r.t. the 1st camera

%% Kanatani {R,h}
% [U,~,~] = svd(Essential);
% h2 = U(:,3);
% 
% K = -skew_sym(h2)*Essential;
% 
% [U,~,V] = svd(K);
% 
% R2 = U*diag([1 1 det(U*V')])*V';

% Choose the one rotation and translation that makes sense
vis_pts1 = []; vis_pts2 = [];
for k=1:4
    R1 = eye(3);
    R2 = rot(:,:,k);
    h1 = [0 0 0]';
    h2 = t(:,:,k);
    
    %     figure, hold on;
    %     plot3(h1(1),h1(2),h1(3),'rx')
    %     xlabel('x'), ylabel('y'), zlabel('z')
    %     axis equal, grid on, hold on
    %     plot3(e1(1),e1(2),e1(3),'cx')
    %     plot3(h2(1),h2(2),h2(3),'bx')
    for i=1:size(v_k,2)
        vis_pts1(1:3,i) = v_k(:,i);
        vis_pts2(1:3,i) = R2'*u_k(:,i) - R2'*h2;
        %         plot3(vis_pts1(1,i), vis_pts1(2,i), vis_pts1(3,i),'b.')
        %         plot3(vis_pts2(1,i), vis_pts2(2,i), vis_pts2(3,i),'g.')
    end
    
    plot3(h1(1),h1(2),h1(3),'rx')
    xlabel('x'), ylabel('y'), zlabel('z')
    axis equal, grid on, hold on
%     plot3(e1(1),e1(2),e1(3),'cx')
%     plot3(e2(1),e2(2),e2(3),'gx')
    plot3(h2(1),h2(2),h2(3),'bx')
    for i=1:10:size(v_k,2)
        %         plot3(vis_pts1(1,i), vis_pts1(2,i), vis_pts1(3,i),'b.')
        %         plot3(vis_pts2(1,i), vis_pts2(2,i), vis_pts2(3,i),'g.')
        plot3([h1(1) vis_pts1(1,i)],[h1(2) vis_pts1(2,i)],[h1(3) vis_pts1(3,i)],'b')
        plot3([h2(1) vis_pts2(1,i)],[h2(2) vis_pts2(2,i)],[h2(3) vis_pts2(3,i)],'g')
    end

     %% Triangulate
    r = zeros(3,size(v_k,2));
    for i=1:size(v_k,2)
         P1 = R1*v_k(:,i) + h1; P0 = h1;
        u = P1-P0;
        Q1 = R2*u_k(:,i) + h2; Q0 = h2;
        v = Q1-Q0;
        
        w0 = P0 - Q0;
        
        a = u'*u;
        b = u'*v;
        c = v'*v;
        d = u'*w0;
        e = v'*w0;
        
        sc = (b*e - c*d) / (a*c - b*b);
        tc = (a*e - b*d) / (a*c - b*b);
        
        s1 = P0 + sc*u;
        s2 = Q0 + tc*v;
        
        r(1:3,i) = (s1 + s2)/2;
    end
    
    %% Triangulate 2 - don't use, wrong points estimated (kept here for documentation purposes)
%     r = zeros(3,size(v_k,2));
%     for i=1:size(v_k,2)
%         A2 = R1*v_k(:,i) + h1; A1 = h1;
%         p1 = A2-A1;
%         B2 = R2*u_k(:,i) + h2; B1 = h2;
%         p2 = B2-B1;
%         
%         % Triangulation finding smallest line segment
%         W = cross(p1,p2);
%         s1 = dot(cross(B1-A1,p2),W)/dot(W,W)*p1;
%         s2 = dot(cross(B1-A1,p1),W)/dot(W,W)*p2;
%         
%         r(1:3,i) = (s1 + s2)/2;
%     end    
    
    %% Triangulate Hartley - ??? Not the right points (kept here for documentation purposes)
%     r3 = triangulate( C1_h,C2_h,[R2 h2] );

    
     %% Triangulate kanatani - needs error correction
    for i=1:size(v_k,2)
        Z(i) = dot(skew_sym(h2)*R2*C2_h(:,i), skew_sym(C1_h(:,i))*R2*C2_h(:,i)) / ...
               norm(skew_sym(C1_h(:,i))*R2*C2_h(:,i))^2;
        r4(1:3,i) = Z(i)*C1_h(:,i);
    end
    
%     for i=1:100%size(v_k,2)
%         plot3(r(1,i),r(2,i),r(3,i),'r.')
%     end
    
    dummy = [0 0 1]';
    dum1 = R1*dummy + h1;
    dum2 = R2*dummy + h2;
    %% Write to ply------------------------------------------------------------
    fileIDw = fopen('truck.ply','w');
    fprintf(fileIDw, '%s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n',...
        'ply', ...
        'format ascii 1.0', ...
        'element vertex 72 # number of points', ...
        'property float x', ...
        'property float y', ...
        'property float z', ...
        'property uchar diffuse_red', ...
        'property uchar diffuse_green', ...
        'property uchar diffuse_blue', ...
        'end_header');
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h1(1),h1(2),h1(3), 255,0,0);
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h2(1),h2(2),h2(3), 0,255,0);
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', dum1(1),dum1(2),dum1(3), 0,255,255);
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', dum2(1),dum2(2),dum2(3), 0,255,255);
    
    for i=1:size(r,2)
        fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', r(1,i),r(2,i),r(3,i), ...
            image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),1), ...
            image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),2), ...
            image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),3));
    end
    fclose(fileIDw);
    
    pause
    close all
    
end

% figure, hold on,
% plot3(vis_pts1(1,:),vis_pts1(2,:),vis_pts1(3,:),'r*');
% plot3(vis_pts2(1,:),vis_pts2(2,:),vis_pts2(3,:),'b*');

%% Results from Ted's reconstruction---------------------------------------
% R1 = [9.9949596630e-01 -3.1473172934e-02 4.1536402695e-03;
% 3.1238581861e-02 9.9836461439e-01 4.7877424173e-02;
% -5.6537019166e-03 -4.7723538507e-02 9.9884458227e-01];
%
%
% R2 = [9.9469112137e-01 -2.9944940031e-03 -1.0286207306e-01;
% 7.3924122568e-03 9.9907332725e-01 4.2400931856e-02;
% 1.0263978424e-01 -4.2936229304e-02 9.9379150475e-01];
%
% h1 = [2.9926314242e-01 6.6521162716e-01 1.2885279692e+00]';
% h2 = [-1.3557032831e+00 5.5810026181e-01 1.4522007657e+00]';

%% Triangulation-----------------------------------------------------------
% for point of intersection between space lines see Kanatani p.108 eq.4.55
r = zeros(3,size(v_k,2));
for i=1:size(v_k,2)
    A2 = R1*v_k(:,i) + h1; A1 = h1;
    p1 = A2-A1;
    B2 = R2*u_k(:,i) + h2; B1 = h2;
    p2 = B2-B1;
    
    % Triangulation finding smallest line segment
    W = cross(p1,p2);
    s1 = dot(cross(B1-A1,p2),W)/dot(W,W)*p1;
    s2 = dot(cross(B1-A1,p1),W)/dot(W,W)*p2;
    
    r(1:3,i) = (s1 + s2)/2;
end

% figure, hold on, grid on, axis equal % xlim([-2 2]), ylim([-2 3]), zlim([-2 2])
% plot3(r(1,:),r(2,:),r(3,:),'r*')
%
% Let's visualize the 3D points w.r.t. the two cameras!
for i=1:size(C1_h,2)
    plot3(h1(1),h1(2),h1(3),'rx')
    xlabel('x'), ylabel('y'), zlabel('z')
    axis equal, grid on, hold on
    plot3(h2(1),h2(2),h2(3),'bx')
    plot3(r(1,i),r(2,i),r(3,i),'r.')
    plot3([h1(1) vis_pts1(1,i)],[h1(2) vis_pts1(2,i)],[h1(3) vis_pts1(3,i)],'b')
    plot3([h2(1) vis_pts2(1,i)],[h2(2) vis_pts2(2,i)],[h2(3) vis_pts2(3,i)],'g')
    %     pause
    %     close all
end

dummy = [0 0 1]';
dum1 = R1*dummy + h1;
dum2 = R2*dummy + h2;
%% Write to ply------------------------------------------------------------
fileIDw = fopen('truck.ply','w');
fprintf(fileIDw, '%s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n',...
    'ply', ...
    'format ascii 1.0', ...
    'element vertex 918 # number of points', ...
    'property float x', ...
    'property float y', ...
    'property float z', ...
    'property uchar diffuse_red', ...
    'property uchar diffuse_green', ...
    'property uchar diffuse_blue', ...
    'end_header');
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h1(1),h1(2),h1(3), 255,0,0);
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', h2(1),h2(2),h2(3), 0,255,0);
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', dum1(1),dum1(2),dum1(3), 0,255,255);
fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', dum2(1),dum2(2),dum2(3), 0,255,255);

for i=1:size(r,2)
    fprintf(fileIDw,'%8.4f %8.4f %8.4f %d %d %d\n', r(1,i),r(2,i),r(3,i), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),1), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),2), ...
        image1(floor(C1_hom_pix(2,i)),floor(C1_hom_pix(1,i)),3));
end
fclose(fileIDw);


