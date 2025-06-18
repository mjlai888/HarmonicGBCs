function T = mkcc(V,T);
%        T = mkcc(V,T)
% If necessary, this function re-orders in counter clockwise order
% the vertices of each triangle in the given triangulation 
    m = size(T,1);
    for i = 1:m
      if triarea(V(T(i,1),:),V(T(i,2),:),V(T(i,3),:)) < 0
        T(i,:) = [T(i,1),T(i,3),T(i,2)];
      end;
    end; 