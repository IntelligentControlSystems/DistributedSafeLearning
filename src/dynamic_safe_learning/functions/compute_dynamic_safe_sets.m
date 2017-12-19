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

function safe_sets = compute_dynamic_safe_sets(system, x0, tau)
%COMPUTE_DYNAMIC_SAFE_SETS computes the dynamic safe sets for the passed
% system. The computations are the one described in the article "Safe 
% Learning for Distributed Systems with Bounded Uncertanties". The input is
% a structure with the following fields
% INPUT:  A,     state matrix
%         B,     input-to-state matrix
%         G,     noise-to-state matrix
%         Q,     ellipsoidal noise shape matrix
%         q,     ellipsoidal noise level-set
%         n,     agents state dimension (column vector)
%         m,     agents input dimension (column vector)
%         p,     agents noise dimension (column vector)
%         H,h,   state constraints Hx <= h  
%         O,o,   input constraints Ox <= o
%         nh,    number of state constraints per system
%         no,    number of input constraints per system
%         ---------------------------------------------
%         x0,    initial condition
% The ouput structure containts
% OUTPUT: P,     safe sets
%         L,     safe control laws
%         Gamma, dynamic update for safe sets
%         alpha, base level-set for safe sets


% extracting the variables form the structure
[~,~,~,~,q,H,h,O,o,n,~,p,~,~,neighbours] = system_2_variables(system);

% system dimension 
M = length(n);

% creating structures with all the local variables
[Ai, Bi, Qi, Gi, ~] = global_2_local_system(system);

% yalmip variables
alpha   = sdpvar(M,1);

xi      = sdpvar(M,1);
for i = 1:M
    Ei{i}   = sdpvar(n(i) , n(i));
    Ki{i}   = sdpvar(size(Bi{i},2) , size(Ai{i},2), 'full');
    Phii{i} = sdpvar(size(Ai{i},2) , size(Ai{i},2));
end
for i = 1:M
    ENi{i} = blkdiag(Ei{neighbours{i}});
end

% transformation for the phii variables 
cumulated_system_length = [0; cumsum(n)];
for i = 1:M
    Tlphi = zeros(sum(system.n), size(Ai{i},2));
    cumulated_local_system_length = [0; cumsum(system.n(system.neighbours{i}))];
    for j = 1:length(neighbours{i})
        Tlphi(cumulated_system_length(j)+1:cumulated_system_length(j+1),...
              cumulated_local_system_length(j)+1:cumulated_local_system_length(j+1)) = eye(system.n(i));
    end
    Phig{i}= Tlphi * Phii{i} * Tlphi';
end

% extracting initial conditon
for i = 1:M
    x_0{i} = x0(cumulated_system_length(i)+1:cumulated_system_length(i+1));
end

% optimization problem constraints
constraints = [];

constraints = [constraints; sum(alpha) == 1; ...
               sum(cat(3,Phig{:}),3) == zeros(cumulated_system_length(end))];

for i = 1:M
    
    initial_condition = [alpha(i) x_0{i}'; x_0{i}  Ei{i}] >= 0;
    
    S_procedure = [(Phii{i} + tau(i)*ENi{i}) , zeros(size(Ai{i},2),p(i)) , (ENi{i}*Ai{i}' + Ki{i}'*Bi{i}'); ...
                    zeros(p(i),size(Ai{i},2)) , xi(i) * Qi{i} , Gi{i}';...
                    Ai{i}*ENi{i} + Bi{i}*Ki{i} , Gi{i} , Ei{i}] >= 0;
    
    constraints = [constraints; ...
                   alpha(i) >= 0; Ei{i} >= 1e-9; ...
                   alpha(i) - tau(i) - xi(i) * q(i) >= 0; ...
                   tau(i) >= 0; xi(i) >= 0; ...
                   initial_condition, S_procedure];
               
    for j = 1:size(H{i},1)
        constraints = [constraints;...
                          [h{i}(j)^2, H{i}(j,:)*ENi{i};...
                           ENi{i}*H{i}(j,:)', ENi{i}] >= 0
                      ];
    end
    
    for j = 1:size(O{i},1)
        constraints = [constraints;...
            [o{i}(j)^2, O{i}(j,:)*Ki{i};...
            Ki{i}'*O{i}(j,:)', ENi{i}] >= 0
            ];
    end
end

% optimization problem: objective
objective = 0;
for i = 1:M
    objective = objective - logdet(Ei{i});
end

% optimization: solution
options = sdpsettings('solver', 'sdpt3');
diagnostics = optimize(constraints, objective, options);

if diagnostics.problem ~= 0
    error('Problem infeasible');
end

% extrapolating safe sets, feedback etc. from optimization variables
for i = 1:M
    safe_sets.Pi{i}     = inv(value(Ei{i}));
    safe_sets.Li{i}     = value(Ki{i})*inv(value(ENi{i}));
    safe_sets.Gammai{i} = inv(value(ENi{i}))*value(Phii{i})*inv(value(ENi{i}));
    safe_sets.alphai{i} = value(alpha(i));
end

end

