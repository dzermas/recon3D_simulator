% Kanatani's methodology to compute motion parameters from the
% Fundamental/Essential matrix and the correspondences

function [R, h, error] = motion_parameters(G, v_c, u_c)

% error: the average residual of the epipolar constraint for all points
% used to compute the fundamental/essential matrix
% Find motion parameters [h,R]
% Kanatani p.340
[U,~,~] = svd(G*G');
h = U(:,3);

sum = 0;
ss = 0;
for i=1:length(v_c)
    ss = ss + u_c(:,i)'*G*v_c(:,i);
    %assert (u_c(:,i)'*G*v_c(:,i) < 10e-03, 'Error Fundamental');
    sum = sum + triple(h, v_c(:,i), G*u_c(:,i));
end
error = ss / length(v_c);
% determine sign of h
if (sum < 0) h = -h; end

K = skew_sym(-h)*G;
[Uk,~,Vk] = svd(K);

R = Uk * diag([1 1 det(Uk*Vk')]) * Vk';