function [cams] = readBundler(file)

fid = fopen(file);

% skip the first line
fgets(fid);

fscanf(fid, '%d', 2);

cams = [readCamera(fid) readCamera(fid) readCamera(fid)];

fclose(fid);

end

function [cam] = readCamera(fid)

f = fscanf(fid, '%f', 1);
r = fscanf(fid, '%f', [2 1]);
R = fscanf(fid, '%f', [3 3])';
t = fscanf(fid, '%f', [3 1]);

cam = struct( ...
    'focal', f, ...
    'radial', r, ...
    'rotation', R, ...
    'center', -R'*t);
    %'center', R * [1 1 1]' + t);


end