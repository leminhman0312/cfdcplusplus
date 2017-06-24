function [p]=laplace2D(p,y,dx,dy,l1norm_target)
    l1norm=1;
    j=size(p,2)
    
while l1norm>l1norm_target
    pn=p;   
    for i=2:(size(p,1)-1)
        p(i,j-1) = (dy^2*(pn(i+1,j)+pn(i,j))+dx^2*(pn(i,j+1)+pn(i,j)))/(2*(dx^2+dy^2));
        p(:,1)=0;
        p(1,:)=p(2,:);
        p(size(p,2),:)=p(size(p,1)-1,:);
        end
l1norm = (sum(abs(p(:,j))) - sum(abs(p(:,j-1)))) /  sum(abs(p(:,j)))
j=j-1
end