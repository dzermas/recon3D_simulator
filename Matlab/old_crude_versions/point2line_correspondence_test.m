%% Find points that correspond to lines
% For image 1
l=1;
for i=1:size(lines1,1)
    % find the points closest to the i_th line in image1
    res = abs(pts1(2,:) - p1(i,1)*pts1(1,:) - p1(i,2));
    ind = find(res < 20);    
    if length(ind) > 1 % at least 3 points on te line
        % Project points onto line
        p_n = p1(i,:)'/norm(p1(i,:),2);
        Pn = eye(2) - p_n*p_n';
        proj_pts = Pn*pts1(1:2,ind) + p1(i,2);
        % Check which points belong to the line segment (p3 = tp1+(1-t)p2 , 0<t<1)
        keep = [];
        for k=1:size(proj_pts,2)
            t = norm(proj_pts(:,k) - [lines1(i,3); lines1(i,4)],2)/ ...
                 norm([lines1(i,1); lines1(i,2)] - [lines1(i,3); lines1(i,4)],2);
            if 0 < t && t < 1
                keep = [keep; ind(k)];
            end
        end
        % If points in line segent exist, create a structure to save all the info
        if ~isempty(keep)
            inlineseg1(l).p = p1(i,:);
            inlineseg1(l).line = i;
            inlineseg1(l).pts_ind = keep;
            inlineseg1(l).pts = pts1(:,keep);
            inlineseg1(l).n_pts = size(pts1(:,keep),2);
            l = l+1;
        end
    end
end

% For image 2
l=1;
for i=1:size(lines2,1)
    % find the points closest to the i_th line in image1
    res = abs(pts2(2,:) - p2(i,1)*pts2(1,:) - p2(i,2));
    ind = find(res < 8);    
    if length(ind) > 1 % at least 3 points on te line
        % Project points onto line
        p_n = p2(i,:)'/norm(p2(i,:),2);
        Pn = eye(2) - p_n*p_n';
        proj_pts = Pn*pts2(1:2,ind) + p2(i,2);
        % Check which points belong to the line segment (p3 = tp1+(1-t)p2 , 0<t<1)
        keep = [];
        for k=1:size(proj_pts,2)
            t = norm(proj_pts(:,k) - [lines2(i,3); lines2(i,4)],2)/ ...
                 norm([lines2(i,1); lines2(i,2)] - [lines2(i,3); lines2(i,4)],2);
            if 0 < t && t < 1
                keep = [keep; ind(k)];
            end
        end
        % If points in line segent exist, create a structure to save all the info
        if ~isempty(keep)
            inlineseg2(l).p = p2(i,:);
            inlineseg2(l).line = i;
            inlineseg2(l).pts_ind = keep;
            inlineseg2(l).pts = pts2(:,keep);
            inlineseg2(l).n_pts = size(pts2(:,keep),2);
            l = l+1;
        end
    end
end
