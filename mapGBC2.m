function [xxf,yyf]=mapGBC2(V,myGBC)
xx=myGBC{1}; 
xxf=zeros(size(xx));
yyf=zeros(size(xx));
for i=1:size(V,1)
    xxf=xxf+myGBC{i}*V(i,1);
    yyf=yyf+myGBC{i}*V(i,2);
end
end