function W=BdryMap(P)
%This is a function taking an input of a polygon P with boundary vertices in P
%and return a set of boundary points  of a square S. 
xmin=min(P(:,1)); xmax=max(P(:,1));
ymin=min(P(:,2)); ymax=max(P(:,2));
M=size(P,1); m1=round(M/4)+1;
W=[];
x=[xmin:(xmax-xmin)/(m1):xmax]; y=ymin*ones(size(x)); x=x(:);y=y(:);
W=[W;x,y];
y=[ymin:(ymax-ymin)/(m1):ymax]; x=xmax*ones(size(y)); x=x(:);y=y(:);
W=[W;x y];
x=[xmin:(xmax-xmin)/(m1):xmax]; y=ymax*ones(size(x)); x=x(:);y=y(:);
W=[W;x(end:-1:1),y];
y1=[ymin:(ymax-ymin)/(m1):ymax]; x1=xmin*ones(size(y1)); x1=x1(:);y1=y1(:);
W=[W;x1, y1(end:-1:1)]; %W=unique(W,'rows');
m=size(W,1); [m, M]
if m>M 
    I=randperm(M); J=I(1:m-M); W(J,:)=[];
end
m=size(W,1); [m,M]

