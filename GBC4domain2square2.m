function [FV1,FV2,color,S2,V,T]=GBC4domain2square2(NofI,Level,W2)
%This  code is written by Dr. Ming-Jun Lai based on his codes for harmonic 
%GBCs and Dr. Tsung-Wei Hu's codes for image deformation and boundary. 
%%Read image in grey scale
I=imread(NofI); %I=imread('Normalcase1.jpg');
I = rgb2gray(I);
% Get level set of high signal, say 150
C=Level; %C=160;
[~,Vertices,Objects]=isocontour(I,C);
%Vertices: all used points. Objects: each cell save a line segment
%We only want the outter line, find a point in highest y, and get the line
%include that poins.
[~,id]=min(Vertices(:,1));
A=cellfun(@(x) ismember(id,x),Objects);
j=find(A);

figure, subplot(121),  imshow(I), hold on;
Points=Objects{j};
plot(Vertices(Points,2),Vertices(Points,1),'Linewidth',2,'Color','r');hold on;
axis([1 512, 1 512])

P=Vertices(Points,:); P1=P(1:10:end,:); P2=[P1(:,2),512-P1(:,1)];
[m,I]=min(P2(:,1)); S2=[P2(I:end,:);P2(1:I-1,:)];
n=size(S2,1); B=[1:n]; B=B'; C=[]; H={};d=10;eps=1e-5;
[V,T]=Triangulation2D(S2,B,H,C,d,eps); T=mkcc(V,T);
subplot(122), triplot(T,V(:,1),V(:,2)), axis image, axis([ 1 512, 1 512])
[xx,yy]=meshgrid(1:512,1:512);
Ind=insideVT(V,T,xx,yy);
nxx=nan*ones(512,512);nyy=nan*ones(512,512);
I=find(Ind==1);
nxx(I)=xx(I);nyy(I)=yy(I);
W=BdryMap(S2); 
myGBC1=FindGBC(S2,nxx(:),nyy(:));
[xx1,yy1]=mapGBC2(W,myGBC1);
FV1=[xx1,yy1]; %mapping to a rectangular domain.

P=imread(NofI); %P=imread('Normalcase1.jpg');
A=double(P(512:-1:1,:,1));color(:,1)=A(:)/255;
A=double(P(512:-1:1,:,2));color(:,2)=A(:)/255;
A=double(P(512:-1:1,:,3));color(:,3)=A(:)/255;
figure, subplot(121), image(P), colormap gray(256)
subplot(122), scatter(nxx(:),nyy(:),5,color), axis([1 512 1 512])
figure, subplot(121), scatter(nxx(:),nyy(:),5,color), axis tight
subplot(122), scatter(yy1,xx1,5,color), axis tight

d=10; n=size(T,1); P=[]; 
for i=1:n
    v1=V(T(i,1),:);v2=V(T(i,2),:); v3=V(T(i,3),:);
    [x,y]=domain_pts(v1,v2,v3,d);
    P=[P;x,y];
end
SV3=unique(P,'rows'); size(SV3)
%[nx,ny,S2]=centernormalize2(SV3(:,1),SV3(:,2),S2); SV3=[nx(:),ny(:)];
myGBC2=FindGBC(S2,SV3(:,1),SV3(:,2));
[xx2,yy2]=mapGBC2(W2,myGBC2); %mapping into a square domain W. 
FV1=[xx2,yy2]; 
FV2=SV3;