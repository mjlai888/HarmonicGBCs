function Mat=buildL(d)
%This matlab computes a matrix of entries of Bernstein polynomials over domain 
%points. 
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
I = [0:d]';
J = d - I;
m = (d+1);
IM = diag(I)*ones(m,m);
JM = diag(J)*ones(m,m);
Mat = (IM/d).^(IM').*(JM/d).^(JM');
IF = gamma(I+1);
JF = gamma(J+1);
A = factorial(d)*ones(m,m)*diag(1./(IF.*JF));
Mat = A.*Mat;
