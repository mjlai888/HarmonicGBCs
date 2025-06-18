function [FV1,W,color]=GBC4domain2square(ImageName,Level)
%This  code is written by Dr. Ming-Jun Lai based on his codes for harmonic 
%GBCs and Dr. Tsung-Wei Hu's codes for image deformation and boundary extraction. 
%%Read image in grey scale
I=imread(ImageName);
I = rgb2gray(I);
C=Level; %C=160;
[~,Vertices,Objects]=isocontour(I,C);
%Vertices: all used points. Objects: each cell save a line segment
%We only want the outter line, find a point in highest y, and get the line
%include that poins.
[~,id]=min(Vertices(:,1));
A=cellfun(@(x) ismember(id,x),Objects);j=find(A);

figure, subplot(121),  imshow(I), hold on;
Points=Objects{j};
plot(Vertices(Points,2),Vertices(Points,1),'Linewidth',2,'Color','r');hold on;
P=Vertices(Points,:); n=size(P,1); 

if n>=1500
    P1=P(1:10:end,:);
else if n>1000
P1=P(1:6:end,:); 
else if n>500
        P1=P(1:3:end,:);
else if n>300
    P1=P(1:2:end,:);
else if n>150
        P1=P;
end
end
end
end
end
P2=[P1(:,2),512-P1(:,1)];
[m,I]=min(P2(:,1)); S2=[P2(I:end,:);P2(1:I-1,:)];
P2=S2;
%P=Vertices(Points,:); P1=P(1:10:end,:); P2=[P1(:,2),512-P1(:,1)];
n=size(P2,1); B=[1:n]; B=B'; C=[]; H={};d=10;eps=1e-5;
[V,T]=Triangulation2D(P2,B,H,C,d,eps); T=mkcc(V,T);
subplot(122), triplot(T,V(:,1),V(:,2)), axis image, axis([ 1 512, 1 512])
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
figure, subplot(121), image(P), colormap gray(256)
subplot(122), scatter(nxx(:),nyy(:),5,color), axis([1 512 1 512])
figure, subplot(121), scatter(nxx(:),nyy(:),5,color), axis tight
subplot(122), scatter(yy1,xx1,5,color), axis tight
FV1=[yy1,xx1]; %mapping to a rectangular domain.
%FVold=[nxx,nyy]; [size(nxx),size(nyy)]