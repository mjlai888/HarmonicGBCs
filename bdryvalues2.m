function g=bdryvalues2(xx,yy,V,i)
%This code is written by Dr. Ming-Jun Lai on Dec. 30, 2018.
k=size(V,1); C=mean(V);
V=[V(k,:);V;V(1,:)];
i=i+1;
v1=V(i-1,:);v2=V(i,:);v3=V(i+1,:);
g=zeros(size(xx));
[b1,b2,b3]=bary(C,v2,v1,xx,yy);
%I=find(abs(b1)<=1e-5);
I=find(abs(b1)<=1e-5&b2>=0 &b3>=0);
g(I)=b2(I);
[b1,b2,b3]=bary(C,v3,v2,xx,yy);
%I=find(abs(b1)<=1e-5);
I=find(abs(b1)<=1e-5&b2>=0 &b3>=0);
g(I)=b3(I);
