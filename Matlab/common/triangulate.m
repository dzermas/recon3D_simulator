function [ XP ] = triangulate( pts1,pts2,P2 )

n=size(pts1,2);
X=zeros(4,n);
for i=1:n
    A=[-1,0,pts1(1,i),0;
        0,-1,pts1(2,i),0;
        pts2(1,i)*P2(3,:)-P2(1,:);
        pts2(2,i)*P2(3,:)-P2(2,:)];
  [~,~,va] = svd(A);
  X(:,i) = va(:,4);
end
XP(:,:,1) = [X(1,:)./X(4,:);X(2,:)./X(4,:);X(3,:)./X(4,:); X(4,:)./X(4,:)];

end