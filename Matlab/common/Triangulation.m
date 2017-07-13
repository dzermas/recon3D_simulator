function X = Triangulation(P1, P2, u, v)

%% Triangulation
% Go to each correspondence and compute the 3D point X (3xN) matrix
X = zeros(3,size(u,2));
for i = 1 : size(u,2)
    % Construct A matrix
    A = [Vec2Skew([u(:,i); 1]) * P1;
         Vec2Skew([v(:,i); 1]) * P2];
    % Solve linear least squares to get 3D point
    [~,~,V] = svd(A'*A);
    point_3d = V(:,end);
    X(:,i) = point_3d(1:3) ./ norm(point_3d(4));
end