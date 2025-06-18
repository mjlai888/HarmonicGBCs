function DerBcoeff = dirder(Bcoeff,lam1,lam2,lam3);
%        DerBcoeff = dirder(Bcoeff,lam1,lam2,lam3);
% Bcoeff = the Bnet coeffs over some triangle of a polynomial of 
% degree d. DerBcoeff = the Bnet coeffs of the directional derivative
% in the direction given by the vector having tcoords [lam1,lam2,lam3].
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
m = size(Bcoeff,1);
d = degree(m);
DerBcoeff = d*de_cast_step(lam1,lam2,lam3,Bcoeff);

function Bout = de_cast_step(lam1,lam2,lam3,Bin);
%        Bout = de_cast_step(lam1,lam2,lam3,Bin);
% This function does one step of the de Casteljau algorithm
% Note: This function assumes that the input arguments are all column vectors
m = size(Bin,1);
d = degree(m);
n = length(lam1);
[I,J,K] = indices(d);
[I1,J1,K1] = indices(d-1);
indx1 = locate(I1+1,J1,K1,I,J,K);
indx2 = locate(I1,J1+1,K1,I,J,K);
indx3 = locate(I1,J1,K1+1,I,J,K);
if size(Bin,2) == 1
   Bout = Bin(indx1)*lam1' + Bin(indx2)*lam2' + Bin(indx3)*lam3';
else
   Bout = Bin(indx1,:)*spdiags(lam1,0,n,n) + Bin(indx2,:)*spdiags(lam2,0,n,n) + ...
      Bin(indx3,:)*spdiags(lam3,0,n,n);
end;


function d = degree(m);
%        d = degree(m)
% d = the degree of the polynomial corresponding to the vector 
% with length m
d = (-3 + sqrt(8*m+1))/2;

function Index = locate(I1,J1,K1,I,J,K);
%        Index = locate(I1,J1,K1,I,J,K);
% Index = the location of the index [I1,J1,K1] in the list 
% [I,J,K]
Index = ismember([I,J,K],[I1,J1,K1],'rows');
Index = find(Index);

