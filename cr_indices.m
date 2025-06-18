function [I1,I2] = cr_indices(r,d);
%        [I1,I2] = cr_indices(r,d);
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
I1 = [];
Start = 1;
D = d + 1;
for j = 0:r
   for k = 0:(r-j)
      new_col = [(Start+k):(Start+k+d-r)]';
      I1 = [I1,new_col];
   end;
   Start = Start + D;
   D = D-1;
end;
I2 = [(-r*r/2+r*(d+3/2)+1):(-r*r/2+r*(d+1/2)+d+1)]';

   