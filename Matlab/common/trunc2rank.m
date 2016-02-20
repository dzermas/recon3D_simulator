function MT=trunc2rank(M,k) 
%%trunc Matrix M to rank(k)
[u,d,v] = svd(M);
sz =size(d); 
d(k+1:sz(1),:) =0 ; 
MT=u*d*v'; 