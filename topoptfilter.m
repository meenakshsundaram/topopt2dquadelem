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
function [filtvar]=topoptfilter(topvar,isnondesign,totele,neigbors)
%
%Input:
%   topvar          - the actual variables for the optimization
%   isnondesign     - an array consisting of ones and zeros such that if
%                     it has a zero it indicates presence of design
%                     domain, else it indicates presence of non-design
%                     domain
%   totele          - the total number of elements
%   neighbors       - the structure containing details of the neigbors
                       % 1. numneig    - number of Neighbors
                       % 2. e          - element numbers of them
                       % 3. distances  - filter radius - Distance to them
                       % 4. divisor    - sum of the above term
%Output:
%   filtvar         - the filtered variables that is the physical 
%                     set of variables
%

filtvar=zeros(totele,1);

for ele=1:totele
    if(isnondesign(ele))
        filtvar(ele)=1;
    else
        filtvar(ele)=sum((neigbors(ele).distances).*topvar(neigbors(ele).e))/neigbors(ele).divisor;
    end
end

