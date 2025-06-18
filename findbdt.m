function B = findbdt(T,V,E,TE,EV);
%        B = findbdr(T,V,E,TE,EV);
% This function finds the boundary vertices for the above triangulation
% and lists them so that any two consecutive vertices are adjacent on
% the boundary.
%% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
 BE = find(sum(TE)==1);
 B = [];
%
% BE = a list of bdr edges, but not yet in cc order.
%
while ~isempty(BE)
    enum = BE(1);
    BEE = E(BE,:);
    e = E(enum,:);
    t = find(TE(:,enum));
    Tri = T(t,:);
    V1 = V(Tri(1),:);
    V2 = V(Tri(2),:);
    V3 = V(Tri(3),:);
    P = (V1 + V2 + V3)/3;
    VL = V(e(1),:);
    VR = V(e(2),:);
    if triarea(VL,VR,P) > 0
       B1 = e';
    else
       B1 = [e(2);e(1)]; 
    end;
    I1 = [1];
    loop = 0;
    i = 1;
    while ~loop
      v1 = B1(i);
      v2 = B1(i+1);
%
% [v1,v2] is the last bdr edge, with the cc direction from v1 to v2 
% 
      next_e = find(((v2==BEE(:,1))&(v1~=BEE(:,2)))|((v2==BEE(:,2))&(v1~=BEE(:,1))));
      if BEE(next_e,1) == v2 
        v3 = BEE(next_e,2);
      else
        v3 = BEE(next_e,1);
      end;
      if v3 ~= B1(1);
         B1 = [B1;v3];
         I1 = [I1;next_e];
         i = i+1;
      else
         I1 = [I1;next_e];
         loop = 1;
      end;
    end;   
    BE = del(BE,BE(I1));
    B = newcol(B,B1);
  end 