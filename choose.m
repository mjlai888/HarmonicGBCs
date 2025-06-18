function y = choose(n,m)
%        y = choose(n,m)
% This function returns the number of ways of choosing n elements from a
% set of m elements. This function assumes that both n and m are nonnegative
% integers.
y = prod(1:m)/(prod(1:n)*prod(1:(m-n)));