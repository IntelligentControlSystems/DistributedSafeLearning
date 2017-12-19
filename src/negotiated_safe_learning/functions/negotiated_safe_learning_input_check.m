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

function negotiated_output = negotiated_safe_learning_input_check(system, safe_sets, x, u)
%NEGOTIATED_SAFE_LEARNING_INPUT_CHECK checks if the provided input is safe to
% apply to the system or not.
% INPUT:  system,    structure with the system,
%         safe_sets, structure with the safe sets
%         x,         current state (with neighbours)
%         u,         input to apply
%         i,         system need to be checked
% OUTPUT: negotiated_output, structure containing
%           delta,   amout of learning input used
%           alphat,  set size
%           input,   safe input

% extracting the variables form the structure
[~,~,~,~,~,~,~,~,~,n,m,~,~,~,neighbours] = system_2_variables(system);

% system dimension 
M = length(n);

% creating structures with all the local variables
[Ai, Bi, Qi, Gi, qi] = global_2_local_system(system);

% extracting safe sets variables from structure
[Pi, Li, ~, alphai] = safe_sets_2_variables(safe_sets);

% extract state with neighbours
cumulated_system_length = [0; cumsum(n)];
for i = 1:M
    T = zeros(sum(system.n), size(Ai{i},2));
    cumulated_local_system_length = [0; cumsum(system.n(system.neighbours{i}))];
    for j = 1:length(neighbours{i})
        T(cumulated_system_length(j)+1:cumulated_system_length(j+1),...
              cumulated_local_system_length(j)+1:cumulated_local_system_length(j+1)) = eye(system.n(i));
    end    
    xNi{i} = T * x;
end

% extract input from global vector
cumulated_input_length = [0; cumsum(m)];
for i = 1:M
    ui{i} = u(cumulated_input_length(j)+1:cumulated_input_length(j+1));
end

% yalmip variables
sigma  = sdpvar(M,1);
delta  = sdpvar(M,1);
alphat = sdpvar(M,1);
for i = 1:M
    S{i} = sdpvar(sum(n(neighbours{i})),sum(n(neighbours{i})));
end

% constraints
constraints = [sum(cat(3,S{:}),3) == 0];
for i = 1:M
    
    % build input g(x,u)
    % SELECT THE RIGHT INPUT AMONG ALL OF THEM!!!
    g = delta(i) * ui{i} + (1 - delta(i)) * Li{i} * xNi{i};     
    
    % local constraints
    constraints = [constraints; 0 <= delta(i) <= 1; ...
                   sigma(i) >= 0; alphat(i) >= 0;...
                   alphai{i} + xNi{i}'*S{i}*xNi{i} >= alphat(i);...
                  [sigma(i) * Qi{i} , zeros(size(Qi{i},1),1), Gi{i}';...
                   zeros(1,size(Qi{i},1)), alphat(i) - sigma(i)*qi{i}, (Ai{i}*xNi{i} + Bi{i}*g)';...
                   Gi{i},(Ai{i}*xNi{i} + Bi{i}*g), inv(Pi{i}) ] >=0];       
end
               
% objective function
objective = -sum(delta);

% solving the feasbility problem
% optimization: solution
options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
diagnostics = optimize(constraints, objective, options);

% if error block the optimization
if diagnostics.problem
    error('Something went wrong in the opimization!');
end

% storing variables
for i = 1:M
    negotiated_output.delta{i}  = value(delta(i));
    negotiated_output.alphat{i} = value(alphat(i));    
    negotiated_output.input{i}  = value(delta(i) * ui{i} + (1 - delta(i)) * Li{i} * xNi{i});   
end
    
end

