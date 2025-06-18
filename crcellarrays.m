function [I1,I2] = CrCellArrays(d,r);
%        [I1,I2] = CrCellArrays(d,r);
% This function returns the cell arrays of indices needed to implement the C^r smoothness 
% across an edge. Note that both I1 and I2 will be 3 X 1 cell arrays.
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
[J1,J2] = cr_indices(r,d);
D1 = zeros(d+1,1);
D2 = zeros(d+1,1);
s1 = d+1;
s2 = 1;
for j = 1:(d+1)
   D1(j) = s1;
   s1 = s1 + d+1-j;
   D2(j) = s2;
   s2 = s2 + d + 2 - j;
end;
I2{1} = flipud(J2);
Temp = D1 - r;
I2{2} = flipud(Temp(1:(d+1-r)));
Temp = D2 + r;
I2{3} = Temp(1:(d+1-r));
I1{1} = J1;
Temp = zeros(size(J1));
for j = 0:r
   Temp(:,j+1)=D1((j+1):(d+1-r+j));
end;
loc = r+2;
back = r+1;
for j = 1:r
   for k = 0:(r-j)
      Temp(:,loc) = Temp(:,loc-back)-1;
      loc = loc+1;
   end;
   back = back-1;
end;
I1{2}=Temp;
Temp = zeros(size(J1));
for j = 0:r
   Temp(:,j+1)=D2((j+1):(d+1-r+j));
end;
loc = r+2;
back = r+1;
for j = 1:r
   for k = 0:(r-j)
      Temp(:,loc) = Temp(:,loc-back)+1;
      loc = loc+1;
   end;
   back = back-1;
end;
I1{3} = flipud(Temp);



