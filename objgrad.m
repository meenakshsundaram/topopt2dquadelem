%
%    This file is a part of the topopt2dquadelem a matlab library. This 
%    is free software: you can redistribute it and/or modify it under 
%    the terms of the GNU Lesser General Public License as published by the 
%    Free Software Foundation, either version 3 of the License, or 
%    (at your option) any later version.
%
%    topopt2dquadelem is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with topopt2dquadelem.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2010 Meenakshi Sundaram   
%
function [Obj,Grad]=objgrad(totele,nodcoord,nodconnect,...
                               totnondesignvec,fixeddofs,force,Dmat,topvar)

%  Computes the objective and gradient for total strain energy
%   Inputs:
%   totele             - Total Number of Elements
%   nodcoord           - Nodal Coordinate
%   nodconnect         - Nodal Connectivity
%   totnondesignvec    - Elements in non design domain
%   fixeddofs          - Fixeddofs
%   force              - Force Vector
%   Dmat               - D Matrix
%   topvar             - Topology Optimization 
%                              Variable (Spurious Density Variable)
    
%   Outputs:
%   Obj                - The Objective Function
%   Grad               - The Gradient of the Objective Function
    
    % -----------------Objective Function--------------------------------------

    %Element Degrees of Freedom
    eledof=zeros(8,1);
    
    %Total number of nodes
    totnod=size(nodcoord,1);
    
    %Total degrees of freedom
    totdofs=2*totnod;

    %The indices for sparse matrix assembly
    X=zeros(totele*64,1);Y=zeros(totele*64,1);Z=zeros(totele*64,1);
    
    %Triplet counter
    ntriplets=0;
    
    %Assembly
    %For every element 
    for ele = 1: totele
        
        %find the nodal connectivity
        nodcon=nodconnect(ele,:);
        
        %find the elemental degrees of freedom
        eledof([1 3 5 7])=2*nodcon-1;
        eledof([2 4 6 8])=2*nodcon;
    
        %The nodal coordinates
        XX=nodcoord(nodcon,:);
        
        %The interpolation of the selection field
        xval=topvar(ele)^3;
        
        % ELement Stiffness Matrix
        Ke=stiffness(Dmat*xval,XX);
        
        %Assembling in row, column, value format
        for i = 1 : 8
            for j = 1:8
    
                ntriplets=ntriplets+1;
                X(ntriplets)=eledof(i);
                Y(ntriplets)=eledof(j);
                Z(ntriplets)=Ke(i,j);
                
            end
        end
        
    end
    
    %Form the sparse matrix
    K=sparse(X,Y,Z,totdofs,totdofs);
    
    %Boundary Conditions
    alldofs=1:totdofs;
    freedofs=setdiff(alldofs,fixeddofs);

    %Solve for displacements
    U=zeros(totdofs,1);
    U(freedofs)=K(freedofs,freedofs)\force(freedofs);
    
    %Calculate the Strain Energy
    Obj=U'*K*U;

    %----------------------Gradient--------------------------------------------
    
    %Gradients
    %The actual set of elements taking part in the design
    actele=setdiff(1:totele,totnondesignvec);

    %The number of actual elements
    nactele=length(actele);

    %Initialise the gradients
    Grad=zeros(nactele,1);

    %Counter to increment and fill up the gradient values
    incr=0;
    
    for ele = actele
        
        %increment the counter
        incr=incr+1;
        
        %the nodal connectivity of the element
        nodcon=nodconnect(ele,:);
        
        %Element degrees of freedom
        eledof([1 3 5 7])=2*nodcon-1;
        eledof([2 4 6 8])=2*nodcon;
    
        %Finding the nodal coordinates
        XX=nodcoord(nodcon,:);
        
        %Differentiation of the interpolated variable
        dxval=3*topvar(ele)^2;
        
        %Derivative of Element Stiffness Matrix
        dKe=stiffness(Dmat*dxval,XX);
    
        %Derivative w.r.t topvar corresponding to that element
        Grad(incr)=-U(eledof)'*dKe*U(eledof);    
        
    end    