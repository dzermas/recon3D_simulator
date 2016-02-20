%% Enforce SIFT points to be computed on the ends of line segments
% Import images
image1 = imread('../../data/trucks/truck1.jpg');
image2 = imread('../../data/trucks/truck2.jpg');

% Grayscale for SIFT features
I1 = single(rgb2gray(image1));
I2 = single(rgb2gray(image2));

%% Call the C code for the LSD!
lines1 = call_LSD('truck1');
lines2 = call_LSD('truck2');

%% Compute SIFT features near the line segment ends!
fc1 = [lines1(:,1)';lines1(:,2)';ones(1,length(lines1))*10;zeros(1,length(lines1))];
[F1,D1] = vl_sift(I1,'frames',fc1,'orientations');

fc2 = [lines2(:,1)';lines2(:,2)';ones(1,length(lines2))*10;zeros(1,length(lines2))];
[F2,D2] = vl_sift(I2,'frames',fc2,'orientations');

% Find correspondences
[MATCHES,SCORES] = vl_ubcmatch(D1, D2);

pts1 = [F1(1,MATCHES(1,:)); F1(2,MATCHES(1,:)); ones(1, length(MATCHES))];
pts2 = [F2(1,MATCHES(2,:)); F2(2,MATCHES(2,:)); ones(1, length(MATCHES))];

% Plot images and their feature points
figure, imshow(image1), hold on
plot(F1(1,MATCHES(1,:)),F1(2,MATCHES(1,:)),'r*'), hold off

figure, imshow(image2), hold on
plot(F2(1,MATCHES(2,:)),F2(2,MATCHES(2,:)),'b*'), hold off

for i=1:10:length(MATCHES)
    figure(99), imshow([image1 image2]), hold on
    plot(F(1,MATCHES(1,i)),F(2,MATCHES(1,i)),'r*'),
    plot((size(image1,2) + F2(1,MATCHES(2,i))),F2(2,MATCHES(2,i)),'b*'),
    plot([F(1,MATCHES(1,i)) (size(image1,2) + F2(1,MATCHES(2,i)))], [F(2,MATCHES(1,i)) F2(2,MATCHES(2,i))])
%     keep = [keep input(prompt)];
    pause, close 99
end