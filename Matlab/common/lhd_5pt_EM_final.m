%%last modify: 31/12/2005: Prepare for ICPR Paper:
%%hongdong LI
%%_ICPR is a good version, which
%% corrects many error's in old lhd_5pt_EM_matlab code!!!

%%re-organized version for ICPR-2006.
%%23-Nov-2005 by lhd
%%remove all evalf(), which may cause numerical problem.
%%WHY three methods get all threee different results??
% use two views and five point to solve for EM
%coded by hongdong li, 2005.
% clearvars;
% close all;
% format long;
% digits(4) ;

function [EMatrices] = lhd_5pt_EM_final(matches) 

EMatrices = [];  

format long;
digits(40) ;

%number of matches desired
% no_matches = 5;
% foc = 1.00;
% noise_std_percentage = 0;

% %%generate 3D data
% X = rand(3,no_matches);
% for i = 1:no_matches
%     X(1,i) = (X(1,i) * 4*foc -2*foc ) ;
%     X(2,i) = (X(2,i) * 4*foc -2*foc );
%     X(3,i) =  X(3,i)*foc + foc;
% end
% tm = ones(no_matches,1);
% true_X = [X', tm ] ; %% hongdong li , output true 3D points
% % X = X ;
% %% scatter3(X(1,:),X(2,:),X(3,:));
% true_X(:,4)  = tm;
% 
% %%%%%%%%%% generate translation and rotation
% t = 2*rand(3,1)*foc;
% t = t./norm(t); %%%% make the norm(EEt) =1.0 for testing purpose
% true_t = t;
% T = [0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0];
% 
% 
% rotation_multplier = 30;
% theta = 2/360 * 2 * pi * rand * rotation_multplier;
% n = 2/360 * 2 * pi * rand * rotation_multplier;
% p = 1/360 * 2 * pi * rand * rotation_multplier;
% R(1,1) = (1 - cos(p)) * cos(n)* ( cos(n) * cos(theta) + sin(theta) * sin(n) ) + cos(p)* cos(theta);
% R(1,2) = (1 - cos(p))* cos(n) * ( sin(n) *cos(theta) - sin(theta) *cos(n) ) - cos(p) *sin(theta);
% R(1,3) = sin(n) *sin(p);
% R(2,1) = (1 - cos(p)) *sin(n)  *( cos(n) *cos(theta) + sin(theta)* sin(n) ) + cos(p)* sin(theta);
% R(2,2) = (1 - cos(p)) *sin(n) * ( sin(n) *cos(theta) - sin(theta) *cos(n) ) + cos(p)* cos(theta);
% R(2,3) = -cos(n) * sin(p);
% R(3,1) = -sin(p) * ( sin(n) * cos(theta) - sin(theta) * cos(n));
% R(3,2) = sin(p) * ( cos(n) * cos(theta) + sin(theta) * sin(n));
% R(3,3) = cos(p);
% 
% 
% %%%%%%generate Projection MAtrices
% f1 = foc;
% f2 = foc;
% 
% % P1 = zeros(3,4);
% K1 = diag([f1,f1,1]);
% P1= K1*[eye(3),zeros(3,1)]  ;
% 
% % P2 = zeros(3,4);
% K2 = diag([f2,f2,1]);
% P2= K2*[R,t] ;
% t = t;
% T=T;
% R=R;
% true_E = T*R;
% true_E =true_E/norm(true_E) ;
% 
% true_F = inv(K2)*true_E*inv(K1) ;
% true_F=true_F/norm(true_F) ;
% 
% perfect_matches = zeros(no_matches,4);
% x1 = zeros(no_matches,1);
% y1 = zeros(no_matches,1);
% 
% x2 = zeros(no_matches,1);
% y2 = zeros(no_matches,1);
% 
% %%%%%%%%begin projection onto two images
% img1 = zeros(3,1);
% img2 = zeros(3,1);
% for i = 1:no_matches
%     img1 =  P1*true_X(i,:)' ;
%     img2 =   P2*true_X(i,:)' ;
%     x1(i) = img1(1)/img1(3);
%     y1(i) = img1(2)/img1(3);
%     x2(i) = img2(1)/img2(3);
%     y2(i) = img2(2)/img2(3);
% end
% 
% perfect_matches(:,1) = x1(:);
% perfect_matches(:,2) = y1(:);
% 
% perfect_matches(:,3) = x2(:);
% perfect_matches(:,4) = y2(:);
% 
% im_size = max(x1(:))-min(x1(:))
% 
% matches =perfect_matches;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Begin the real work ...

%% now compute the null space corresponding to the first five  points:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q =  [ matches(1:5,1)' ;matches(1:5,2)'; 1, 1, 1, 1, 1]' ;
qp = [ matches(1:5,3)' ;matches(1:5,4)'; 1, 1, 1, 1, 1]' ;

M  = zeros(5,9); %%get the five points
for i = 1:5
    M (i,:) = [qp(i,1)*q(i,1) , qp(i,1)*q(i,2),qp(i,1)*q(i,3), qp(i,2)*q(i,1), qp(i,2)*q(i,2), qp(i,2)*q(i,3),qp(i,3)*q(i,1) ,qp(i,3)*q(i,2) , qp(i,3)*q(i,3) ] ;
end
% M = M ;
%%now compute the null space:
N= null(M);

%%begin symbolic computation
syms x y z w e E EQU
e = x*N(:,1) +y*N(:,2) +z*N(:,3)+ N(:,4) ;

E = transpose(reshape(e,3,3));

equ(1) = det(E);
EQU =  simplify(2*E*transpose(E)*E - trace(E*transpose(E))*E);
equ(2:10) = EQU ;
equ = transpose(equ);

for i =1:10
    equ(i) = maple('expand',equ(i));
    equ(i) = maple('collect', equ(i),'[x,y]','distributed' );
    equ(i) = maple('sort',equ(i),'[x,y,z]');
end
maple('evalf', equ);

%pause
%% [5,9,12,13 ] = z^3,z^2,z,1
%%others = [x,y]
oth_idx =[1,2,3,4,6,7, 8,10,11] ;
for i =1:10
    for j =1:9 %% the number of operand is 13.
        oper = maple('op', oth_idx(j), equ(i)) ;
        coef = maple('coeffs',oper,'[x,y]');
        CoefMtx(i,j) = coef ;
    end
    %%coeffs for z,
    CoefMtx(i,10) = maple('op', 5, equ(i)) + maple('op', 9, equ(i))+ maple('op', 12, equ(i))+ maple('op', 13, equ(i)) ;
end


disp('Now compute the determinant of the 10x10 coefMtx...');
EQ = det(CoefMtx)
% EQ = maple('evalf', maple('sort', EQ))
disp('now solving the Det equation...');
EQ= maple('evalf', EQ);
aaa= solve(EQ);


%%%now compute x,y
real_root_count =0 ;

estimated_EM = zeros(3,3,10);

for i=1:10  %%test all ten solutions of z
    if isreal(aaa(i)),
        z = aaa(i) ;
        cc =subs(CoefMtx) ;
        cc=trunc2rank(cc,9);
        xy_root = null(cc);
        xy_root = xy_root./xy_root(10);
        real_root_count  =real_root_count +1 ;
        x= double(xy_root(8));
        %x = sqrt(xy_root(5))
        %x = double(xy_root(1))^(1/3)
        y= double(xy_root(9));
        %y = sqrt(xy_root(7))
        %y = double(xy_root(4))^(1/3)
        
        e = x*N(:,1) +y*N(:,2) +z*N(:,3)+ N(:,4) ;
        est_E = double(reshape(e,3,3)');
        est_E =est_E./norm(est_E) ;
        estimated_EM(:,:,real_root_count)  = est_E ;
        err_in_EM_est = norm( abs(est_E) - abs(true_E))/norm(true_E);
        
    end
end

number_of_real_root = real_root_count
Groundtruth_EM = true_E
disp('**** all the estimated real Essential Matrices (up to a +/- sign) *******');

for i =1:real_root_count
    EM = estimated_EM(:,:,i)
end



%%note: Hartley normalization is essential, unless im_size == 2.0 .

