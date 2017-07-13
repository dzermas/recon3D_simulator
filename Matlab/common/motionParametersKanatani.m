% Find motion parameters [R,h]
% Kanatani p.340
% Kanatani's methodology to compute motion parameters from the
% Fundamental/Essential matrix and the correspondences
function [R, h] = motionParametersKanatani(G, v_c, u_c)
% INPUT:
% G: Fundamental/Essential matrix
% v_c : 3xN homogeneous correspondences of viewpoint 1
% u_c : 3xN homogeneous correspondences of viewpoint 2
% OUTPUT
% R, h : rotation and translation

[U,~,~] = svd(G*G');
h = U(:,3);

% Assuming depths (Z, Z') from cameras 1 and 2 are positive,
% determine sign of h
sum = 0;
for i=1:length(v_c)
    sum = sum + sign(triple(h, v_c(:,i), G*u_c(:,i)));
end

if (sum < 0)
    h = -h;
end

K = skewSym(-h)*G;
[Uk,~,Vk] = svd(K);

R = Uk * diag([1 1 det(Uk*Vk')]) * Vk';

% Now actually see if depths are positive
Zsum = 0;
for i=1:length(v_c)
    Z = dot(cross(h,R*u_c(:,i)), cross(v_c(:,i),R*u_c(:,i)));
    Zsum = Zsum + Z;
end

if Zsum > 0
    return
else
    h = -h;
end
