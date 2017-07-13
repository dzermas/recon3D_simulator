function F_new = correct_F (C1_h, C2_h, F)

%% It doesn't work as expected!

% Assuming that the measurements we get are correct, let's optimize for F
X = zeros(9, size(C1_h,2));
for i=1:size(C1_h,2)
    X(:,i) = [C2_h(1,i)*C1_h(1,i)
              C2_h(1,i)*C1_h(2,i)
              C2_h(1,i)*C1_h(3,i)
              C2_h(2,i)*C1_h(1,i)
              C2_h(2,i)*C1_h(2,i)
              C2_h(2,i)*C1_h(3,i)
              C2_h(3,i)*C1_h(1,i)
              C2_h(3,i)*C1_h(2,i)
              C2_h(3,i)*C1_h(3,i)];
      
    error(i) = C2_h(:,i)' * F * C1_h(:,i);
end

Df = X' \ error';
DF = reshape(Df, 3,3)';

F_new = F - DF;

% for i=1:size(C1_h,2)
%     
%       
%     error_newnew(i) = C2_h(:,i)' * FF_new * C1_h(:,i);
% end

[U,S,V] = svd(F_new);
F_new = U*diag([S(1,1) S(2,2) 0])*V';