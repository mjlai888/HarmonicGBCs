function Kg = stiffM(V,T,d);
%    K=stiffM(V,T,d);
%This code is written by Clayton Mersmann with a help and codes from Dr.
%Ming-Jun Lai in Oct. 2017.
%If you find these programs useful for your research, please
%acknowledge Dr. Ming-Jun Lai's generocity by citing his paper
%"Awanou, G., Lai, M. J. and Wenston, P., 
%The Multivariate Spline Method for
%Scattered Data Fitting and Numerical Solution 
%of Partial Differential Equations, in Wavelets and Splines, 
%Nashboro Press, (2006) 
%edited by G. Chen and  M. J. Lai, pp. 24--74."

%        could pass areavector av from mass matrix if we wanted to
n = size(T,1); m=(d+2)*(d+1)/2;
Vl=V(T(:,:)',:);

idx3=[1:1:n]'*3; idx2=idx3-1; idx1=idx2-1;

E1=Vl(idx1,:)-Vl(idx3,:); E2=Vl(idx2,:)-Vl(idx3,:); 
E1c=[E1,zeros(size(idx1))]; E2c=[E2,zeros(size(idx1))];
ecr=cross(E1c,E2c);
Av=sqrt(sum(ecr.*ecr,2))/2;

ridx1=meshgrid(1:3,1:2);
cidx1=meshgrid(1:2,1:3)'; 
rvidx=repmat(ridx1(:),n,1);
cvidx=repmat(cidx1(:),n,1);
addrow=3*floor([0:1:6*n-1]/6)';
addcol=2*floor([0:1:6*n-1]/6)'; 

Bm1=sparse(rvidx+addrow,cvidx+addcol,Vl')';
sumbary=[1,1,1]; Idn=speye(n,n);
B=sparse([Bm1;sparse(kron(Idn,sumbary))]);
vxv=repmat([1;0],n,1); vyv=repmat([0;1],n,1);
sv=ones(n,1); vz=zeros(2*n,1);

bvdx=B\[vxv;sv];
bvdy=B\[vyv;sv]; 
av=B\[vz;sv]; %directional coordinates for origin

thdx=bvdx-av; 
thdy=bvdy-av; 

Mat1=build(d-1); %for helmholtz will also need Mat=build(d)
B=indices2d(d); C=indices2d(d-1);
I=C(:,1);J=C(:,2);K=C(:,3);
Index1=find(ismember(B,[I+1,J,K],'rows')); %same as find(ismember(B,[I+1,J,K,L],'rows')
Index2=find(ismember(B,[I,J+1,K],'rows'));
Index3=find(ismember(B,[I,J,K+1],'rows'));

Bin=eye(m);
S1=Bin(Index1,:); S2=Bin(Index2,:); S3=Bin(Index3,:);

%Global Dx 
Dxg=d*sparse(kron(spdiags(thdx(idx1),0,size(idx1,1),size(idx1,1)),S1)+...
    kron(spdiags(thdx(idx2),0,size(idx2,1),size(idx2,1)),S2)+...
    kron(spdiags(thdx(idx3),0,size(idx3,1),size(idx3,1)),S3)); 
%Dxg=d*sparse(kron(sparse(diag(thdx(idx1))),S1)+kron(sparse(diag(thdx(idx2))),S2)+...
   % kron(sparse(diag(thdx(idx3))),S3));  
%Global Dy
Dyg=d*sparse(kron(spdiags(thdy(idx1),0,size(idx1,1),size(idx1,1)),S1)+...
    kron(spdiags(thdy(idx2),0,size(idx2,1),size(idx2,1)),S2)+...
    kron(spdiags(thdy(idx3),0,size(idx3,1),size(idx3,1)),S3)); 

%Dyg=d*sparse(kron(sparse(diag(thdy(idx1))),S1)+kron(sparse(diag(thdy(idx2))),S2)+...
%    kron(sparse(diag(thdy(idx3))),S3));

Mat1gV=sparse(kron(spdiags(Av,0,size(Av,1),size(Av,1)),Mat1));
Kg=(Dxg'*Mat1gV*Dxg + Dyg'*Mat1gV*Dyg);

