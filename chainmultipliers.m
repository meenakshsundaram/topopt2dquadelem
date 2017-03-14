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
function [ChainMul]=chainmultipliers(neigbors,totele,isnondesign)
%  To find the derivative of the physical variable with respect to the
%  actual variable
%  Input:
%   neighbors       - the structure containing details of the neigbors
                       % 1. numneig    - number of Neighbors
                       % 2. e          - element numbers of them
                       % 3. distances  - filter radius - Distance to them
                       % 4. divisor    - sum of the above term
%   totele          - total number of elements
%   isnondesign     - vector of the size of the number of elements
%                     indicating whether 1 the element should be designed
%                     or not 0 
%  Output:
%   ChainMul        - derivative multipliers 
%  
    %Initialization of the Chain_Mul array
    ChainMul=zeros(totele);
    
    %For every element put in the value
    for ele=1:totele
        if(~isnondesign(ele))
            ChainMul(ele,neigbors(ele).e)=neigbors(ele).distances/(neigbors(ele).divisor);
        end
    end
    
    %Obtain only the actual elements
    %Actually vectors
    actele=find(~(isnondesign));
    ChainMul=ChainMul(actele,actele);