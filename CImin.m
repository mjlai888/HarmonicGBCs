function [c,GBs] = CImin(A,B,F,G,epsilon,max_it,tol)
%        c = CImin(A,B,F,G,epsilon,max_it,tol)
% CONSTRAINED ITERATIVE MINIMIZATION:
% This function uses iteration to minimize 1/2*c'*A*c-c'*F
% subject to the constraint B*c = G.
%If you find these programs useful for your research, please
%acknowledge Dr. Ming-Jun Lai's generocity by citing his paper
%"Awanou, G., Lai, M. J. and Wenston, P., 
%The Multivariate Spline Method for
%Scattered Data Fitting and Numerical Solution 
%of Partial Differential Equations, in Wavelets and Splines, 
%Nashboro Press, (2006) 
%edited by G. Chen and  M. J. Lai, pp. 24--74."
warning off
Mat = A + 1/epsilon*B'*B;
ds=whos('Mat');
GBs=ds.bytes/1e9; %size of Mat in GB
c = Mat\(F + 1/epsilon*B'*G);
cvg = 0;
it_count = 0;
while ~cvg & it_count <= max_it
    cold = c;
    c = Mat\(A*c + 1/epsilon*B'*G);
    it_count = it_count + 1;
    norm(c - cold,inf);
    cvg = norm(c - cold,inf) <= tol*norm(c,inf);
end
%if ~cvg
    disp('Warning: Constrained iterative minimization failed to converge');
%else
%    disp(['Constrained iterative minimization converged after ',int2str(it_count),' iterations']);
%end;
