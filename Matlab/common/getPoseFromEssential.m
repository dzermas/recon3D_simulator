% This is the "standard" found in Zisserman Hartley's book
function [R, t] = getPoseFromEssential(E, u, v)

% Translation
[U,~,V] = svd(E);

h{1} =  V(:,3);
h{2} = -V(:,3);
h{3} = h{1};
h{4} = h{2};

% Rotation
w = [0 -1 0; 1 0 0; 0 0 1];
z = [0 1 0; -1 0 0; 0 0 1];

rot{1} = U * w * V';
rot{2} = U * z * V';
rot{3} = rot{2};
rot{4} = rot{1};

% Chirailty
current_best = 0;
correct_pose = 0;
for i=1:4
    pos_counter = 0;
    % Build camera matrix 1 and 2
    % P1
    P1 = [eye(3) zeros(3,1)];
    % P2
    P2 = [rot{i} -rot{i}*h{i}];
    X = Triangulation(P1, P2, u, v);
    for j=1:size(X,2)
        if X(3,j) > 0
            pos_counter = pos_counter + 1;
        end
    end
    if pos_counter > current_best
        current_best = pos_counter;
        correct_pose = i;
    end
end

assert(correct_pose ~= 0);

R = rot{correct_pose};
t = h{correct_pose};
    