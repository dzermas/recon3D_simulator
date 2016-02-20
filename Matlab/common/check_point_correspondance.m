%% Early version of the 2 images reconstruction. Contains crude implementations. Kept for reference reasons.

%% Tool that visualizes the correspondences between two images
A = [image1 image2];

figure, imshow(A), hold on

plot(pts1(:,1), pts1(:,2), 'r*')

[x,y] = ginput;

% [neighbors distances] = kNearestNeighbors( XYZsurvey(:,1:2), [ xpick ypick], k );
[neighbors distances] = kNearestNeighbors( pts1(:,1:2), [ x y ], 1);

pts = zeros(length(neighbors),4);
for k = 1:length(neighbors)
    x2 = size(image1,2) + pts2(neighbors(k),1);
    y2 = pts2(neighbors(k),2);
       xy = [x(k) y(k); x2 y2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

% Plot the feature points of the first image
plot(pts1(:,1), pts1(:,2), 'wo')
