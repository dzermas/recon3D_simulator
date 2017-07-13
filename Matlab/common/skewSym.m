function X = skewSym(x)

%% Change a vector into its skew symmetric form. Vector can be column or row.

if (size(x,2)==3)
    x = x';
end

X=[ 0  -x(3)  x(2) ; 
   x(3)  0   -x(1) ; 
  -x(2)  x(1)  0  ];
