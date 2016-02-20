function r = distPointToLineSegment( xy0, xy1, xyP )
% r = distPointToLineSegment( xy0, xy1, xyP )

if( nargin < 3 )
  selfTest();
  r=0;
else
vx = xy0(1)-xyP(1);
vy = xy0(2)-xyP(2);
ux = xy1(1)-xy0(1);
uy = xy1(2)-xy0(2);
lenSqr= (ux*ux+uy*uy);
detP= -vx*ux + -vy*uy;

if( detP < 0 )
  r = norm(xy0-xyP,2);
elseif( detP > lenSqr )
  r = norm(xy1-xyP,2);
else
  r = abs(ux*vy-uy*vx)/sqrt(lenSqr);
end
end


 function selfTest()
   %#ok<*NASGU>
 disp(['invalid args, distPointToLineSegment running (recursive)  self-test...']);

ptA = [1;1]; ptB = [-1;-1];
ptC = [1/2;1/2];  % on the line
ptD = [-2;-1.5];  % too far from line segment
ptE = [1/2;0];    % should be same as perpendicular distance to line

distCtoAB = distPointToLineSegment(ptA,ptB,ptC)
distDtoAB = distPointToLineSegment(ptA,ptB,ptD)
distEtoAB = distPointToLineSegment(ptA,ptB,ptE)
figure(1); clf; 
plot([ptA(1) ptB(1)],[ptA(2) ptB(2)],'r-x'); hold on;
plot(ptC(1),ptC(2),'b+'); plot(ptD(1),ptD(2),'g+'); plot(ptE(1),ptE(2),'k+'); hold off;
axis([-3 3 -3 3]);
end

end