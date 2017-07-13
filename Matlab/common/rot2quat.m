function q = rot2quat( R )
% Transforms a rotation matrix into a quaternion
trace = R(1,1) + R(2,2) + R(3,3);
if (trace > 0)
    s = 0.5 / sqrt( trace + 1 );
    q(1) = 0.25 / s;
    q(2) = ( R(3,1) - R(2,3) ) * s;
    q(3) = ( R(1,3) - R(3,1) ) * s;
    q(4) = ( R(2,1) - R(1,2) ) * s;
elseif ( R(1,1) > R(2,2) && R(1,1) > R(3,3) )
    s = 2 * sqrt( 1 + R(1,1) - R(2,2) - R(3,3) );
    q(1) = ( R(3,2) - R(2,3) ) / s;
    q(2) = 0.25 * s;
    q(3) = ( R(1,2) + R(2,1) ) / s;
    q(4) = ( R(1,3) + R(3,1) ) / s;
elseif ( R(2,2) > R(3,3) )
    s = 2 * sqrt( 1 + R(2,2) - R(1,1) - R(3,3));
    q(1) = ( R(1,3) - R(3,1) ) / s;
    q(2) = ( R(1,2) + R(2,1) ) / s;
    q(3) = 0.25 * s;
    q(4) = ( R(2,3) + R(3,2) ) / s;
else 
    s = 2 * sqrt( 1 + R(3,3) - R(1,1) - R(2,2) );
    q(1) = ( R(2,1) - R(1,2) ) / s;
    q(2) = ( R(1,3) + R(3,1) ) / s;
    q(3) = ( R(2,3) + R(3,2) ) / s;
    q(4) = 0.25 * s;
end

