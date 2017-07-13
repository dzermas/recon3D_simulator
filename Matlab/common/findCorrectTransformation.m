function [R, t] = findCorrectTransformation(rot, t, C1_h, C2_h)

% Input: rot  : 4 rotations
%        t    : 4 translations
%        C1_h : homogeneous pixel coordinates Image 1
%        C2_h : homogeneous pixel coordinates Image 2

R1 = eye(3);
t1 = [0 0 0]';
    
R2 = rot(:,:,1);
t2 = t(:,1);

% Triangulate kanatani
Z.Z1 = [];
for i=1:size(points,2)
	% No Kanatani correction reconstruction
	Z1 = dot(cross(h2,R2*C2_h(:,i)), cross(C1_h(:,i),R2*C2_h(:,i))) / ...
        norm(cross(C1_h(:,i),R2*C2_h(:,i)))^2;
end