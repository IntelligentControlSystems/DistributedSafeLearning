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

function [A,B,G,Q,q,H,h,O,o,n,m,p,nh,no,neighbours] = system_2_variables(system)
%SYSTEM_2_VARIABLES Transforms the system structure in single variables

A         = system.A;
B         = system.B;
G         = system.G;
Q         = system.Q;
q         = system.q;
H         = system.H;
h         = system.h; 
O         = system.O;
o         = system.o;
n         = system.n;
m         = system.m;
p         = system.p;
nh        = system.nh;
no        = system.no;
neighbours = system.neighbours;

end

