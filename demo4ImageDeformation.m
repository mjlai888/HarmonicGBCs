%This code is to deform one medical image to another. 
% It is written by Dr. Ming-Jun Lai based on Dr. Tsung-wei Hu's codes to 
% extract the boundary of the domain of interest. 
%The deformation is achieved by using the bivariate spline implementation
%of GBC functions by using the codes developed by Dr. Tsung-Wei Hu under the
% supervision of Dr. Ming-Jun Lai. The mathematical details can be found in Tsung-Wei Hu's
%dissertation, 2023 and the paper by Hu, Deng, and Lai, 2022. 
%Dr. Ming-Jun Lai rewrote the following for medical image deformation 
%on May, 8, 2024 and revised on June 12, 2025.
%First, we find GBC map based on Image 1
ImageName='Malignant3.jpg'; 
[A,B,C]=GBC4domain2square(ImageName,150);
%Then we find GBC map based on Image 2.
[m,~]=size(B); 
if m>150
    I=randperm(m); J=I(1:150); I=sort(J);
    B1=B(I,:); %be sure to include the four corners.
else    
    B1=B;
end
[A2,B2,C2,S,V,T]=GBC4domain2square2('Normalcase1.jpg',150,B1);
k=knnsearch(A,[A2(:,2),A2(:,1)]);
newcolor=C(k,:);
I=imread('Normalcase1.jpg'); J=imread(ImageName);
figure, subplot(131), image(I), 
subplot(132), image(I), hold on, scatter(B2(:,1),512-B2(:,2),5,newcolor);
subplot(133), image(J);