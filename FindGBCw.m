function myGBC=FindGBCw(V,xx,yy,caseNum) 
%This code is written by Mr. Tsung-wei Hu under Dr. Ming-Jun Lai's
%supervision in 2022. 
n=size(V,1);
Vold=V; m=size(V,1); H={};C=[]; B=[1:n]; B=B'; d=0.1; tol=1e-8;
%[V,T]=Triangulation2D(V,B,H,C,d,tol);
[V,T]=Tdomain2N(Vold,20); %Another way to find a triangulation.
d=8;r=0; maxit=5; eps=1e-4;  %d=3r+2

H=smoothnessM2(V,T,d,r); 
[A,B,C,D,B1,B2,C1] = bnetw(V,T,d,'PDEweights',caseNum);
K=stiffw2P(V,T,d,A,B,C,D,B1,B2,C1); % the weighted stiff matrix
n=size(K,1); Zeros=zeros(n,1);
n=size(H,1); Z4H=zeros(n,1);
[E,TE,~,~,~] = tdataM(V,T); 
myGBC={}; %This saves the function values for each GBC function. 
parfor i=1:m
  [G,B] = dirchletGBC(V,T,TE,E,d,'bdryvalues2',Vold,i);
  c=CImin(K,[H;B],Zeros,[Z4H;G],eps,maxit,1e-6);
  I=find(c<-1e-5);size(I), % number of negative values
  [norm(H*c,inf),norm(B*c-G,inf)] % check the smoothness and boundary constraints
    Z = seval(V,T,c,xx,yy);
    Z(Z<0)=0;
    zz=reshape(Z,size(xx)); myGBC{i}=zz;
end

end
