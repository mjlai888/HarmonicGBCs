function [Edges,TE,TV,EV,Bdr] = tdataM(V,T);
%Altered from tdata(V,T) by Clay Mersmann, Fall 2017
%        [E,TE,TV,EV,B] = tdata(V,T)
% Input = Vertices and triangles of a triangulation.
% Output:
% Edges = a list of edges
% TE = Num(triangles) X Num(edges) matrix whose (i,j)th entry 
% = 1 if edge #j is an edge for triangle #i.
% TV = Num(triangles) X Num(vertices) matrix whose (i,j)th entry
% = 1 if vertex #j is a vertex for triangle #i.
% EV = Num(edges) X Num(vertices) matrix whose (i,j)th entry
% = 1 if edge #i has vertex #j as an endpoint.
% Bdr = a list of the bdr vertice in cc order on the boundary
% BdrEdges = a list of the bdr edges in cc order on the boundary
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

 szt= size(T,1); 
 Tt=T';
 
 E1=sort(T(:,1:2),2); E2=sort(T(:,2:3),2); E3=sort([T(:,3),T(:,1)],2);
 ix3=[1:1:szt]*3; ix2=ix3-1; ix1=ix2-1;
 E=nan(szt*3,2);
 E(ix1,:)=E1; E(ix2,:)=E2; E(ix3,:)=E3;
 
 [Edges,ia,ic]=unique(E,'rows','stable');
 
 tidx=floor(1/3*[0:1:szt*3-1]')+1;
 TE=sparse(tidx,ic,1);
 TV = sparse(tidx,Tt(:),1);
 
 sze=size(Edges,1);
 eidx=floor(1/2*[0:1:sze*2-1]')+1; 
 Et=Edges';
 EV=sparse(eidx,Et(:),1);

 Bdr = findbdt(T,V,Edges,TE);

    
     
      
 
