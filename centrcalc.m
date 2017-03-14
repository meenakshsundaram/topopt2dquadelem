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
function [neigbors,centroids]=centrcalc(nodcoord,nodcon,isnondesign,...
                                       totnondesignvec,totele,filterradius)
%   To find the neighboring set of elements for every element also to 
%   find a certain set of parameters to aid the construction of the 
%   density filter
%
%   Input
%       nodcoord        - Nodal Coordinate
%       nodcon          - Nodal Connectivity
%       isnondesign     - An array consisting of ones and zeros. 
%                         If it is one then it is non-design 
%                         element else it is not
%       totnondesignvec - Set of all elements which are a part of the 
%                         nondesign domain
%       totele          - Total number of elements
%       filterradius    - Radius of the filter within which the neighboring 
%                         elements need to be found
%   Output 
%       neighbors       - the structure containing details of the neigbors
                           % 1. numneig    - number of Neighbors
                           % 2. e          - element numbers of them
                           % 3. distances  - filter radius - Distance to them
                           % 4. divisor    - sum of the above term
%       centroids       - the value of centroid for each element    

    

    %Centroids Computation
    centroids=zeros(totele,2);
    for ele=1:totele
        centroids(ele,:)=sum(nodcoord(nodcon(ele,:),:))/4;
    end

    %Initialize the structure
    neigbors=repmat(struct('numneig',{},'e',{},'distances',{},'divisor',{}),totele,1);

    %For every element compute the structure
    for ele=1:totele
    
        %verify whether it is a part of the non_design element sets
        if(~isnondesign(ele))
        
            %if it doesn't
            cloc=centroids(ele,:);
            
            %find all those guys within the filter radius
            dist=sqrt((cloc(1)-centroids(:,1)).^2 + (cloc(2)-centroids(:,2)).^2);
            ind=find(dist<filterradius);
            
            %remove the guys within the non_design domain
            ind=setdiff(ind,intersect(ind,totnondesignvec));
            
            %find the distance values
            distset=filterradius-dist(ind);
            
            %find their sum
            divisor=sum(distset);
            
            %organise the structure
            neigbors(ele)=struct('numneig',max(size(ind)),'e',ind,'distances',distset,'divisor',divisor);
            
        end
    end