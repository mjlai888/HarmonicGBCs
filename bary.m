function [lam1,lam2,lam3] = bary(V1,V2,V3,X,Y);
%        [lam1,lam2,lam3] = bary(V1,V2,V3,X,Y);
% This function returns the barycentric coordinates with respect to the 
% triangle [V1,V2,V3] of the points (X,Y);
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
%If you find these programs useful for your research, please
%acknowledge Dr. Ming-Jun Lai's generocity by citing his paper
%"Awanou, G., Lai, M. J. and Wenston, P., 
%The Multivariate Spline Method for
%Scattered Data Fitting and Numerical Solution 
%of Partial Differential Equations, in Wavelets and Splines, 
%Nashboro Press, (2006) 
%edited by G. Chen and  M. J. Lai, pp. 24--74."
One = ones(size(X(:)'));
A = [1 1 1 ;V1(1),V2(1),V3(1);V1(2),V2(2),V3(2)];
lam = A\[One;X(:)';Y(:)'];
lam1 = reshape(lam(1,:),size(X));
lam2 = reshape(lam(2,:),size(X));
lam3 = reshape(lam(3,:),size(X));
