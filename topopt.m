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
function [topvar,filtvar,iter]=topopt(nodcoord,nodcon,nondesign,...
                                  filterradius,fixeddofs,force,volfrac,tol)
% topology optimization for minimization of total strain energy given a
% volume constraint. It uses the density filter and 4 node plane stress 
% quad elements with displacements based fem formulation
%
% Input:
%   nodcoord        - Nodal Coordinates
%                       x y
%                       x y
%                       Assumed to be in ascending order of node numbers
%   nodcon          - Nodal Connectivity
%                       It is written in anti-clockwise order
%                       nod1 nod2 nod3 nod4
%                       Assumed to be in ascending order of element numbers
%   nondesign       - A structure for non design containing both filled and 
%                       empty
%                       .fill  -non design elements that are filled
%                       .empty -non design elements that are empty              
%   filterradius    - radius of the density filtering
%   fixeddofs       - the fixed degrees of freedom
%   force           - the force vector
%   volfrac         - the maximum volume fraction available
%   tol             - the tolerance for convergence in the relative norm 
%                     sense
% Output:	
%   topvar          - the design variables
%   filtvar         - the filtered variables 
%                     the variables with physical meaning
%   iter            - the number of iterations taken


    %Initialization
    %total number of elements
    totele=size(nodcon,1);
    
    %Initialization of the selection field
    %called topvar here
    topvar=volfrac*ones(totele,1);

    %adjusting the non design domains
    topvar(nondesign.fill)=1;
    topvar(nondesign.empty)=0.001;
    
    %total set of nondesign vectors
    totnondesignvec=[nondesign.fill; nondesign.empty];

    %Make an arrays of 1 and 0s 
    %if the element belongs to a non-design domain let it be 1 
    %else let it be 0
    isnondesign=zeros(totele,1);
    for ele=1:totele
        if(ismember(ele,totnondesignvec))
            isnondesign(ele)=1;
        end
    end
    
    %iteration change value is set to a number greater than tolerance
    change=1.+tol;
    
    %the iteration counter
    iter=0;
    
    %the Young's Modulus is preferrably scaled to 1
    E=1.;
    
    %the poissons ratio
    nu=0.3;
    
    %The Constitutive Tensor for plane stress
    DMat=E/(1-nu^2)*[   1.0      nu          0.0 ;
                        nu       1.0         0.0 ;
                        0.0      0.0        (1.0-nu)/2.0];
    %The Constitutive Tensor for plane strain
         %(E/((1.0+nu)*(1.0-2.0*nu)))*[ 1.0-nu  nu      0.0;
         %                              nu      1.0-nu  0.0;
         %                              0.0     0.0     (1.0-2.0*nu)/2.0 ];
    %Calculation of neighbors
    [neigbors,~]=centrcalc(nodcoord,nodcon,isnondesign,totnondesignvec,totele,filterradius);
    % neighbors is a structure that contains the following
    % 1. numneig    - number of Neighbors
    % 2. e          - element numbers of them
    % 3. distances  - filter radius - Distance to them
    % 4. divisor    - sum of the above term
    
    %Calculation of chain multipliers
    %this a matrix of the derivative of the physical variable with respect to
    %actual variable
    [chainmul]=chainmultipliers(neigbors,totele,isnondesign);

    %FILTERING ONCE PRIOR TO THE ITERATIONS 
    [filtvar]=topoptfilter(topvar,isnondesign,totele,neigbors);

    %------------PRECOMPUTATION DONE-------------------------------------------

    %Begin Iteration
    while(change>tol && iter<1000)
        
        %Increment the iter counter
        iter=iter+1;
        
        %Store old Val for comparison
        topvarold=topvar;
        
        %Objgrad
        [Obj,Grad]=objgrad(totele,nodcoord,nodcon,totnondesignvec,fixeddofs,force,DMat,filtvar);
        
        %Grad multiplied with chainmul
        Grad=chainmul'*Grad;
        
        %Optimality criteria update
        [topvar,filtvar]=opcriteria(Grad,topvar,isnondesign,totnondesignvec,totele,neigbors,volfrac);
        
        %Detect chanfge
        change=norm((topvarold-topvar),2)/norm(topvarold,2);
        
        % Print output every 10 iterations
        if(mod(iter,10)==0)
            sprintf([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',Obj) ...
                ' Vol.: ' sprintf('%6.3f',sum(filtvar)/(totele)) ...
                ' ch.: ' sprintf('%6.3f',change)])
        end

    end
    % Print the final iteration
    sprintf([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',Obj) ...
           ' Vol.: ' sprintf('%6.3f',sum(filtvar)/(totele)) ...
            ' ch.: ' sprintf('%6.3f',change)])
