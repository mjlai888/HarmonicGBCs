function Iv = indices2d(d)
%        [I,J,K] = indices(d)
% This function computes the index vectors I,J,K associated with our linear ordering
% for the Bnet coefficients of a polynomial of degree d over an arbitrary triangle.
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
I = zeros(m,1);
J = I;
K = I;
Mark = 1;
for j = d:(-1):0
   I(Mark:(Mark+j)) =[j:(-1):0]';
   J(Mark:(Mark+j)) =[0:j]';
   K(Mark:(Mark+j)) =(d-j)*ones(j+1,1);
   Mark = Mark+j+1;
end;

Iv=[I,J,K];