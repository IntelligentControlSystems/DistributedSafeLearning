% SAFE LEARNING
% Copyright (C) 2017 Andrea Carron
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% State matrix
A = [0.5 0.3 0.4  0;...
    0  0.2  0  0.4;...
    0.6  0  0.7 0.2;...
    0  0.6  0  0.1];

% state dimension
n = [2 ; 2];

% Input to state matrix
B = [0 0;...
    1 0;...
    0 0;...
    0 1];

% input dimension
m = [1 ; 1];

% Noise to state matrix
G = [0 0;...
     1 0;...
     0 0;...
     0 1];

% noise dimension
p = [1 ; 1];

% noise ellipsoidal shape
Q = eye(2);

% noise level-sets
q = [1e-3 ; 1e-3];

% system state constraints
H{1} = [1  0  0  0;...
       -1 0  0  0;...
        0  1  0  0;...
        0 -1  0  0];

h{1} = ones(4,1);

H{2} = [0  0  1  0;...
        0  0 -1  0;...
        0  0  0  1;...
        0  0  0 -1];

h{2} = ones(4,1);

% number of state constraints per system
nh = [4; 4];

% input constraints
O{1} = [ 1;...
        -1];
o{1} = [1; 1];

O{2} = [ 1;...
        -1];

o{2} = [1; 1];

% number of input constraints per system
no = [2; 2];

% neigbouring nodes, structure with the neighbouring nodes
% the self loop must be included. The i-th element is related to node i
neighbours{1} = [1 2];
neighbours{2} = [1 2];

% storing all the variables in a structure
system.A  = A;
system.B  = B;
system.G  = G;
system.Q  = Q;
system.q  = q;
system.H  = H;
system.h  = h;
system.O  = O;
system.o  = o;
system.n  = n;
system.m  = m;
system.p  = p;
system.nh = h;
system.no = no;
system.neighbours = neighbours;

% deleting old variables
clearvars -except system