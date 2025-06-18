function A = newcol(B,c);
%        A = newcol(B,c)
% This function sets A = [B,c] with B or c padded with zeros as needed
  [n,k] = size(B);
  m = length(c);
  if n > m
    A = [B,[c;zeros(n-m,1)]];
  else
    A = [[B;zeros(m-n,k)],c];
  end;