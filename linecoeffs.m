function c = LineCoeffs(Mat,d,g,V1,V2,varargin);
%        c = LineCoeffs(Mat,d,g,V1,V2,varargin);
% This function computes the Bnet coefficients over the Line [V1,V2]
% of the poly g(x) of degree <= d
% Call Mat=buildL(d); to get matrix Mat.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
I = [0:d]';
J = d - I;
X = (J*V1(1) + I*V2(1))/d;
Y = (J*V1(2) + I*V2(2))/d;
b = feval(g,X,Y,varargin{:});
c = Mat\b;
