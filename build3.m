function Mat = build3(d1,d2,d3)
%        Mat = build3(d1,d2,d3)
% This function builds the matrix needed to compute the triple product of polynomials 
% of degrees d1, d2, and d3 over  an arbitrary triangle.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
[I1,J1,K1] = indices(d1);
[I2,J2,K2] = indices(d2);
[I3,J3,K3] = indices(d3);
m1 = (d1+1)*(d1+2)/2;
m2 = (d2+1)*(d2+2)/2;
m3 = (d3+1)*(d3+2)/2;
Mat = zeros(m1*m2,m3);
for j1 = 1:m1
   for j2 = 1:m2
       for j3 = 1:m3
           F1 = factorial(d1)/(factorial(I1(j1))*factorial(J1(j1))*factorial(K1(j1)));
           F2 = factorial(d2)/(factorial(I2(j2))*factorial(J2(j2))*factorial(K2(j2)));
           F3 = factorial(d3)/(factorial(I3(j3))*factorial(J3(j3))*factorial(K3(j3)));
           F = factorial(d1+d2+d3)/(factorial(I1(j1)+I2(j2)+I3(j3))*factorial(J1(j1)+J2(j2)+J3(j3))*...
               factorial(K1(j1)+K2(j2)+K3(j3)));
           Mat((j1-1)*m2+j2,j3) = F1*F2*F3/F;
       end;
   end;
end;
Mat = Mat/choose(2,d1+d2+d3+2);
