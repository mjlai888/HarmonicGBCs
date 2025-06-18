function Z = seval(V,T,c,X,Y);
%        Z = seval(V,T,c,X,Y)
% This function evaluates the piecewise polynomials over the triangulation [V,T] defined 
% by the vector c of Bnet coeffs at the points with coords X and Y .
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
k = size(T,1);
Z =nan(size(X));
m = length(c)/k;
tol = 100*eps;
for j = 1:k
   [lam1,lam2,lam3] = bary(V(T(j,1),:),V(T(j,2),:),V(T(j,3),:),X,Y);
   I = find(lam1 >=-tol & lam2 >= -tol & lam3 >= -tol);
   if ~isempty(I)
      Z(I) = loceval(lam1(I),lam2(I),lam3(I),c(((j-1)*m+1):(j*m)));
   end;
end;

function y = loceval(lam1,lam2,lam3,Bcoeff);
%        y = loceval(lam1,lam2,lam3,Bcoeff);
% This function evaluates the poly with Bnet coefficients in the vector Bcoeff 
% at the points with barycentric coordinates lam1, lam2, lam3.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
m = length(Bcoeff);
d = degree(m);
for j = 1:d
   Bcoeff = de_cast_step(lam1,lam2,lam3,Bcoeff);
end;
y = Bcoeff';
   
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
