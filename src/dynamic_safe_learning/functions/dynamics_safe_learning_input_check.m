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

function safe = dynamics_safe_learning_input_check(system, safe_sets, x, u, i)
%DYNAMICS_SAFE_LEARNING_INPUT_CHECK checks if the provided input is safe to
% apply to the system or not.
% INPUT:  system,    structure with the system,
%         safe_sets, structure with the set
%         x,         current state (with neighbours)
%         u,         input to apply
%         i,         system need to be checked
% OUTPUT: safe, 1 if the input is safe, 0 otherwise

% extracting the variables form the structure
[A,B,G,Q,q,H,h,O,o,n,m,p,nh,no,neighbours] = system_2_variables(system);

% system dimension 
M = length(n);

% creating structures with all the local variables
[Ai, Bi, Qi, Gi, qi] = global_2_local_system(system);

% extracting safe sets variables from structure
[Pi, Li, Gammai, alphai] = safe_sets_2_variables(safe_sets);

% yalmip variable
sigma = sdpvar(1,1);

% building the LMI to check and constraint on sigma
constraints = [sigma >=0;...
               [Gi{i}'*Pi{i} *Gi{i} - sigma * Qi{i} ,  Gi{i}'*Pi{i}*(Ai{i}*x +Bi{i}*u);...
               (Ai{i}*x +Bi{i}*u)'*Pi{i}*Gi{i}, (Ai{i}*x +Bi{i}*u)'*Pi{i}*(Ai{i}*x +Bi{i}*u) - alphai{i} - x'*Gammai{i}*x + sigma*qi{i}] <=0];       

% objective function
objective = 0;

% solving the feasbility problem
% optimization: solution
options = sdpsettings('solver', 'mosek', 'verbose', 0);
diagnostics = optimize(constraints, objective, options);

% returing 1 if safe, and 0 otherwise
safe = ~diagnostics.problem;
end

