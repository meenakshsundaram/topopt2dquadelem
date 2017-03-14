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

% Demo script for topology optimization of a bridge type loading
% with a volume constraint of 0.35 of the total volume

%% Initialization
clear;
clc;

% Loading Mesh data
load('Bridge/fixednodes.txt');
load('Bridge/forcenodes.txt');
load('Bridge/nodalconnectivity.txt');
load('Bridge/nodalcoord.txt');
load('Bridge/nondesignelefill.txt');
load('Bridge/nondesigneleempty.txt');
nondesignele=struct('fill',nondesignelefill,'empty',nondesigneleempty);
clear nondesignelefill;
clear nondesigneleempty;
%%

%filter radius
filterradius=10;

%fixed degrees of freedom
fixeddofs=union(2*fixednodes,2*fixednodes-1);

%total number of nodes and elements
totnod=size(nodalcoord,1);
totele=size(nodalconnectivity,1);

%force vector
force=zeros(2*totnod,1);
force(2*forcenodes)=-1.0;

%Maximum volume fraction that can be used
volfrac = 0.35;

%Tolerance
tol=1e-2;

%Calling the topology optimization routine
[topvar,filtvar,iter]=topopt(nodalcoord,nodalconnectivity,nondesignele,filterradius,fixeddofs,force,volfrac,tol);

%% Plotting results
col=[linspace(1,0,101)',linspace(1,0,101)',linspace(1,0,101)'];

h=figure(1);clf;set(gca,'fontsize',32,'linewidth',2);set(h,'Colormap',col);
rho = (filtvar-min(filtvar))/(max(filtvar)-min(filtvar));
col=[linspace(1,0,101)',linspace(1,0,101)',linspace(1,0,101)'];
set(h,'Colormap',col);
hold on
for ele = 1:totele
    nodcon=nodalconnectivity(ele,:);
    cx=nodalcoord(nodcon,1);
    cy=nodalcoord(nodcon,2);
    patch(cx,cy,(1-rho(ele))*[1 1 1],'edgecolor','none');
end
axis equal;axis tight;
plot(nodalcoord(fixednodes,1),nodalcoord(fixednodes,2),'ob','MarkerSize',6','MarkerFaceColor',[0,1,0]);
plot(nodalcoord(forcenodes,1),nodalcoord(forcenodes,2),'ob','MarkerSize',6','MarkerFaceColor',[0,1,0]);


