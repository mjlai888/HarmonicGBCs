function [G,B] = dirchletGBC(V,T,TE,E,d,g,varargin);
%        [G,B] = dirchlet(V,T,TE,E,d,g,varargin);
% This function returns the coeffs in the vector c of the piecewise poly over the 
% triangulation [V,T] corresponding to the Dirchlet boundary conditions defined
% by the function g. 
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
m = (d+1)*(d+2)/2;
Mat=buildL(d);
k = size(T,1);
BE = find(sum(TE)==1);
n = length(BE);
N = n*(d+1);
Indx2 = zeros(N,1);
G = zeros(N,1);
A=[1,2;2,3;3,1;2,1;3,2;1,3];
I1 = crcellarrays(d,0);
for j = 1:n
   LocIndx = [((j-1)*(d+1)+1):(j*(d+1))];
   BT = find(TE(:,BE(j)));
   v1 = E(BE(j),1);
   V1 = V(v1,:);
   v2 = E(BE(j),2);
   V2 = V(v2,:);
   i1 = find(T(BT,:) == v1);
   i2 = find(T(BT,:) == v2);
   e = find(ismember(A,[i1,i2],'rows'));
   if e <=3
      G(LocIndx) = linecoeffs(Mat,d,g,V1,V2,varargin{:});
      Indx2(LocIndx) = (BT-1)*m + I1{e};
   else
      e = e - 3;
      G(LocIndx) = linecoeffs(Mat,d,g,V2,V1,varargin{:});
      Indx2(LocIndx) = (BT-1)*m + I1{e};
   end;
end
[Indx2,I] = unique(Indx2);
G = G(I); 
N = length(Indx2);
ONE = ones(N,1);
Indx1 = [1:N]';
B = sparse(Indx1,Indx2,ONE,N,k*m);
