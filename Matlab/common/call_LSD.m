function [lines, p] = call_LSD(file, threshold)
% lines: x1, y1, x2, y2, line_width, ... see LSD README (we only keep 1,2,3,4 columns)
% p: a, b (as in y=ax+b)

% file: name of the image (e.g. truck1)
% threshold on the length of lines

system(strcat(['./truck_lines.sh ' file]));

fid = fopen(strcat([file '.txt']),'r');
lines_cell = textscan(fid, '%f%f%f%f');
fclose(fid);

lines_init = [lines_cell{1} lines_cell{2} lines_cell{3} lines_cell{4}];
     
% Here we compute the distance of the line segment and make a decision on
% keeping it based on "threshold"
dist = zeros(size(lines_init,1),1);
for i=1:size(lines_init,1)
    dist(i) = sqrt((lines_init(i,2) - lines_init(i,4))^2 + (lines_init(i,1) - lines_init(i,3))^2);
end

lines = lines_init(dist > threshold,:);

% It is vital to know the parameters of the lines as well!
% a.k.a. (a,b) in y=ax+b
p = zeros(size(lines,1),2);
for j=1:size(lines,1)
    p(j,1) = (lines(j,2) - lines(j,4)) / (lines(j,1) - lines(j,3)); % a = (y1-y2)/(x1-x2)
    p(j,2) = lines(j,2) - p(j,1)*lines(j,1); % b = y1-ax1   
end


    