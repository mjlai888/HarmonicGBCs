%This code is to deform one image to another. It is written
%by Dr. Ming-Jun Lai based on Dr. Tsung-wei Hu's codes to extract
%the boundary of the domain of interest. It is written on May, 8, 2024.
%The deformation is achieved by using the bivariate spline implementation
%of GBC functions. The mathematical details can be found in Tsung-Wei Hu's
%dissertation, 2023 and the paper by Hu, Deng, and Lai, 2022.
ImageName='gridlines2.jpg'; 
%The following are two typical domains, L-shaped and H-shaped domains.
L=[0 0;1/2 0;1 0;1 1/2;1/2 1/2;1/2 1;0 1;0 1/2];
L=511*L; L=L+1; %L-shaped domain for images of size 512x512.
[V,T]=Tdomain2N(L,10);
[Edges,TE,TV,EV,B] = tdataM(V,T);
domain=V(B,:);
[L1,L2,B,C]=GBCdomain2square(domain,ImageName);
nxx=L1(:,1:512);nyy=L1(:,513:1024);
figure, subplot(121),scatter(nxx(:),512-nyy(:),5,C), axis([1,512,1,512])
subplot(122), scatter(L2(:,1),L2(:,2),5,C),axis([1,512,1,512])
S2=[0 0;0.5 0;1 0;1 0.5;1 1;1.5 1;2 1;2 0.5;2 0;2.5 0;3 0;3 0.5;3,1;3 1.5;3 2;3 2.5;3 3;2.5 3;...
   2 3;2 2.5;2 2;1.5 2;1 2;1 2.5;1 3;0.5 3;0 3;0 2.5;0 2;0 1.5;0 1;0 0.5];
S2=S2/3*511+1; %H-shaped domain for images of size 512x512.
%Construct harmonic GBC functions of the polygon
[V2,T2] = Tdomain2N(S2,10); 
%figure, triplot(T2,V2(:,1),V2(:,2))
[Edges,TE,TV,EV,B] = tdataM(V2,T2);
domain=V2(B,:);
[H1,H2,B2,C2]=GBCdomain2square(domain,ImageName);
nxx=H1(:,1:512);nyy=H1(:,513:1024);
figure, subplot(121),scatter(nyy(:),512-nxx(:),5,C), axis([1,512,1,512])
subplot(122), scatter(H2(:,1),H2(:,2),5,C2),axis([1,512,1,512])
return
