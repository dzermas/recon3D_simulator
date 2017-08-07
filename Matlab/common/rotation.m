function DR = rotation ( axis, angle )
% rotation matrix around "axis" for "angle" radians

switch axis
    case 'z'
        DR = [cos(angle) -sin(angle) 0;
              sin(angle)  cos(angle) 0;
                  0           0      1];

    case 'y'
        DR = [ cos(angle) 0 sin(angle);
                   0      1     0;
              -sin(angle) 0 cos(angle)];

    case 'x'       
        DR = [ 1     0           0;
               0 cos(angle) -sin(angle);
               0 sin(angle) cos(angle)];
end