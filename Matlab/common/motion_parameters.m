function [R, h] = motion_parameters(G, v_c, u_c)

% Find motion parameters [h,R]
% Kanatani p.340
[U,~,~] = svd(G*G');
h = U(:,3);

sum = 0;
for i=1:length(v_c)
    sum = sum + triple(h, v_c(:,i), G*u_c(:,i));
end

% determine sign of h
if (sum < 0) h = -h; end

K = skew_sym(-h)*G;
[Uk,~,Vk] = svd(K);

W = [0 -1 0;
     1  0 0;
     0  0 1];


R = Uk * diag([1 1 det(Uk*Vk')]) * Vk';