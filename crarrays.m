function [I1,I2] = crarrays(d,r);
%        [I1,I2] = crarrays(d,r);
%If you find these programs useful for your research, please
%acknowledge Dr. Ming-Jun Lai's generocity by citing his paper
%"Awanou, G., Lai, M. J. and Wenston, P., 
%The Multivariate Spline Method for
%Scattered Data Fitting and Numerical Solution 
%of Partial Differential Equations, in Wavelets and Splines, 
%Nashboro Press, (2006) 
%edited by G. Chen and  M. J. Lai, pp. 24--74."for j = 0:r
for j = 0:r
   [J1,J2] = crcellarrays(d,j);
   I1{j+1} = J1;
   I2{j+1} = J2;
end;
end

function [I1,I2] = crcellarrays(d,r);
%        [I1,I2] = CrCellArrays(d,r);
% This function returns the cell arrays of indices needed to implement the C^r
% across an edge condition. Note that both I1 and I2 will be 3 X 1 cell arrays.
[J1,J2] = cr_indices(r,d);
D1 = zeros(d+1,1);
D2 = zeros(d+1,1);
s1 = d+1;
s2 = 1;

%D1 and D2 are indices of other sides of triangle; i.e. not [1:d+1];
for j = 1:(d+1) 
   D1(j) = s1;
   s1 = s1 + d+1-j;
   D2(j) = s2;
   s2 = s2 + d + 2 - j;
end;
I2{1} = flipud(J2); %flipupsidedown
Temp = D1 - r;
I2{2} = flipud(Temp(1:(d+1-r)));
Temp = D2 + r;
I2{3} = Temp(1:(d+1-r));
I1{1} = J1;
Temp = zeros(size(J1));
%I2 cell of local indices in clockwise order?

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
end

function [I1,I2] = cr_indices(r,d);
%        [I1,I2] = Cr_indices(r,d);
% This function returns indices needed to formulate the equations between
% Bnet coefficients that are a consequence of the C^r condition
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
I2 = [(-r*r/2+r*(d+3/2)+1):(-r*r/2+r*(d+1/2)+d+1)]'; %written to obfuscate?
%upper limit is (d+1-r)+...+(d)+(d)...this is the "row" where smoothness 
%condition stops overlapping from triangle to triangle
end

function [I,J,K] = indices(d)
%        [I,J,K] = indices(d)
% This function computes the index vectors I,J,K associated with our linear ordering
% for the Bnet coefficients of a polynomial of degree d over an arbitrary triangle.
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
end