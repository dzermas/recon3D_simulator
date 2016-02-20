function [] = pyramid(P,l,w,h) 
% PYRAMID will accept four inputs. P0 is an array of 3 scalar numbers for 
% the origin (x, y, and z). l is the length of the box in the x-direction, 
% w is the width of the box in the y-direction, and h is the height of the 
% box in the z-direction. The functin will draw a square pyramid. 
% Input: Four inputs, an array of the point of origin, a length, width, and 
% height. 
% Output: A pyramid drawn with a set transparency and different colors for 
% the faces 
% Proecessing: The first thing I will do is use the figure funciton in 
% order to prevent any previous figures from being overwritten.
x = [P(1),P(1)+l,P(1)+l,P(1)]; 
y = [P(2),P(2),P(2)-w,P(2)-w]; 
z = [P(3),P(3),P(3),P(3)]; 
fill3(x, y, z ,'blue'), hold on  
x2 = [P(1),P(1)+l,P(1)+ l/2]; 
y2 = [P(2),P(2),P(2)-w/2]; 
z2 = [P(3),P(3),P(3)+h]; 
fill3(x2, y2, z2,'green'), hold on
x3 = [P(1)+l,P(1)+l,P(1) + l/2]; 
y3 = [P(2), P(2)-w,P(2)- w/2]; 
z3 = [P(3),P(3),P(3)+h]; 
fill3(x3, y3, z3,'red'), hold on 
x4 = [P(1)+l,P(1),P(1)+ l/2]; 
y4 = [P(2)-w,P(2)-w,P(2)- w/2]; 
z4 = [P(3),P(3),P(3)+h];  
fill3(x4,y4,z4,'green'), hold on 
x5 = [P(1),P(1),P(1) + l/2]; 
y5 = [P(2),P(2)-w,P(2)- w/2]; 
z5 = [P(3),P(3),P(3)+h];  
fill3(x5,y5,z5,'red') 