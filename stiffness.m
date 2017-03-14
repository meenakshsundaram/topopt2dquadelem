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
function [K] = stiffness(D,XX)
%  To find the element stiffness matrix for 4 noded plane stress
%  element 
%  Input:
%       D   -   Constitutive Matrix   (3x3)
%       XX  -   Nodal coordinates of the four nodes in counter clockwise
%               direction
%  Output:
%       K   -   Element stiffness matrix
%
 
    %Initializing Element Stiffness Matrices
    K=zeros(8);

    %Thickness of the element
    Thic = 1.0;

    % XG YG are gauss point matrices and WGT is weights
    c=1/sqrt(3);
    XG = [c, c,-c,-c];
    YG = [c,-c, c,-c];
    WGT =[1.0,1.0,1.0,1.0];

    B=zeros(3,8);
    
    % Finding stiffness matrix
    % Loop over Gauss Points And Integrate
    for GSP = 1:4
    
        R = XG(GSP);	
        S = YG(GSP);
        RP = 1.0+R;SP = 1.0+S;RM = 1.0-R;SM = 1.0-S;
    
        P = [0.25*SP -0.25*SP  -0.25*SM  0.25*SM;
             0.25*RP  0.25*RM  -0.25*RM -0.25*RP];
    
         %Jacobian and determinant of the jacobian
        XJ=P*XX;    
        Det = det(XJ);
    
        XJI = [XJ(2,2) -XJ(1,2);
              -XJ(2,1) XJ(1,1)]/Det;
    
        %Strain Displacement Matrix
        B(1,[1 3 5 7]) = XJI(1,:)*P;
        B(2,[2 4 6 8]) = XJI(2,:)*P;
        B(3,[1 3 5 7]) = B(2,[2 4 6 8]);
        B(3,[2 4 6 8]) = B(1,[1 3 5 7]);
        
        %Accumulation of Stiffness Matrix with the weights 
        WT = WGT(GSP)*Thic*Det;
        K=K+WT*B'*D*B;
        
    end
