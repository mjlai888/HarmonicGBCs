function Mat = build(d)
%        Mat = build(d)
% This function builds the matrix needed to compute the inner product of two
% polynomials of degree d over an arbitrary triangle.
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
[I,J,K] = indices(d);
m = (d+1)*(d+2)/2;
Mat = zeros(m,m);
for j = 1:m
   for i = 1:m
      Mat(i,j) = nchoosek(I(j)+I(i),I(j))*nchoosek(J(j)+J(i),J(j))*nchoosek(K(j)+K(i),K(j));
   end;
end;
Mat = Mat/(nchoosek(2*d,d)*nchoosek(2*d+2,2));
