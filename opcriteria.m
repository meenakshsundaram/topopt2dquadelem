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
function [topvarnew,filtvar]=opcriteria(grad,topvar,isnondesign,totnondesignvec,totele,neigbors,volfrac)
%   To perform the Optimality Criteria Update
%   Input:
        % grad              -   gradient of the objective function
        % topvar            -   variables of actual optimization
        % isnondesign       -   an array consisting of ones and zeros. 
        %                       If it is one then it is non-design element 
        %                       else it is not
        % totnondesignvec   - the total non design set of elements
        % totele            - the total number of elements
        % neigbors          - the neigbors structure
        % volfrac           - the max volume fraction it should 
%   Output:
        %topvarnew          - finalised set of actual variables
        %filtvar            - filtered variables (the physical variable set)
    
    %Set the lower and uppser limit    
    l1 = 0; l2 = 100000; 
    
    %Set the move limit and initialize the topvar_new to zeros
    move = 0.2;
    topvarnew=zeros(totele,1);
    
    %copy the non_design variable set
    topvarnew(totnondesignvec)=topvar(totnondesignvec);
    
    %find the actual elements going to take part in the design
    actele=setdiff(1:totele,totnondesignvec);

    %bisection operation
    while (l2-l1 > 1e-4)
        
        %central value
        lmid = 0.5*(l2+l1);
        
        %update
        topvarnew(actele) = max(0.001,...
            max(topvar(actele)-move,min(1.,min(topvar(actele)+move,...
            topvar(actele).*(sign(-grad./lmid).*(abs(-grad./lmid)).^0.4))))...
            );
        
        %filtered variable
        [filtvar]=topoptfilter(topvarnew,isnondesign,totele,neigbors);
  
        %check whether the volfraction has been exceeded and perform the
        %bisection
        if sum(filtvar) - volfrac*totele > 0;
            l1 = lmid;
        else
            l2 = lmid;
        end
    end