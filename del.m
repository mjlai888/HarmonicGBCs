function [lnew,I] = del(l1,l2);
%        [lnew,I] = del(l1,l2)
% This function deletes the set of elements l2 from l1. Note: I = the 
% indices of those elements in l1 which are not in l2
  m = length(l2);
  n = length(l1);
  lnew = l1(:)'; 
  I = 1:n;
  for i = 1:m
    loc = find(l2(i) == lnew);
    if ~isempty(loc)
      k = length(loc);
      for j = 1:k
        J = find(l1 == l2(i));
        I(J) = zeros(size(J));
        lnew = [lnew(1:(loc(j)-1)),lnew((loc(j)+1):length(lnew))];
      end;
    end;
  end;
  I = I(find(I));
  n = size(l1);
  if n(1) > 1
    lnew = lnew';
  end;
     
  
     