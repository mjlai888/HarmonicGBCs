function [X,Y] = domain_pts(V1,V2,V3,d);
[I,J,K] = indices(d);
X = (V1(:,1)*I' + V2(:,1)*J' + V3(:,1)*K')/d;
Y = (V1(:,2)*I' + V2(:,2)*J' + V3(:,2)*K')/d;
X=X(:);
Y=Y(:);
end
