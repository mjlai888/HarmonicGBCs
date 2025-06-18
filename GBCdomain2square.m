function [FV,FV1,W,color]=GBCdomain2square(domain,ImageName)
%This  code is written by Dr. Ming-Jun Lai based on his codes for harmonic 
%GBCs and Dr. Tsung-Wei Hu's codes for image deformation and boundary. 
%%Read image in grey scale
I=imread(ImageName);
I = rgb2gray(I);
[m,n]=size(I(:,:,1));
figure(1), subplot(121),  imshow(I), hold on;
axis([1 m, 1 n])
P1=domain;
P2=[P1(:,2),512-P1(:,1)];
[m,I]=min(P2(:,1)); 
S2=[P2(I:end,:);P2(1:I-1,:)];
P2=S2;
n=size(P2,1); B=[1:n]; B=B'; C=[]; H={};d=10;eps=1e-5;
[V,T]=Triangulation2D(P2,B,H,C,d,eps); T=mkcc(V,T);
figure(1), subplot(122), triplot(T,V(:,2),V(:,1)), axis image, axis([1 512,1 512])
S1=P2;
[xx,yy]=meshgrid(1:512,1:512);
Ind=insideVT(V,T,xx,yy);
nxx=nan*ones(512,512);nyy=nan*ones(512,512);
I=find(Ind==1);
nxx(I)=xx(I);nyy(I)=yy(I);
W=BdryMap(S1); 
myGBC1=FindGBC(S1,nxx(:),nyy(:));
[xx1,yy1]=mapGBC2(W,myGBC1);
P=imread(ImageName); 
A=double(P(512:-1:1,:,1));color(:,1)=A(:)/255;
A=double(P(512:-1:1,:,2));color(:,2)=A(:)/255;
A=double(P(512:-1:1,:,3));color(:,3)=A(:)/255;
figure(2), subplot(121), image(P), colormap gray(256)
subplot(122), scatter(nyy(:),nxx(:),5,color), axis([1 512 1 512])
figure(3), subplot(121), scatter(nyy(:),nxx(:),5,color), axis tight
subplot(122), scatter(xx1,yy1,5,color), axis tight
FV1=[xx1,yy1]; %mapping to a rectangular domain.
FV=[nxx,nyy]; %the original domain