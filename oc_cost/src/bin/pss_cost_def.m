function state_cost = pss_cost_def(x)

%
%   add_state = pss_target_set(x)
%
%   This function allows the user to define any characteristic function
%   in any fashion to define the target set, as long as 
%   the state x (input) can be applied. Return {0,1} based on whether
%   the state belongs to the target set or not.
%   
%       INPUTS: x - state to be considered [x1;x2;...]
%
%       OUTPUT: add_sate - Return {0,1} depending on whether the state 
%			   passed belongs to characteristic function
%
%   Anna Davitian <davitian@ee.ucla.edu>,   CyPhyLab-UCLA 2009

% The user is free to define any characteristic function


state_cost = Inf;


 

 
if (x(1) >= 0 || x(1) <= 2) && ( x(2) >= 0 || x(2) <= 2)
	state_cost = 8;
end
 
if (x(1) > 2 || x(1) <= 4) && ( x(2) > 2 || x(2) <= 4)
	state_cost = 6;
end
 
if (x(1) > 4 || x(1) <= 6) && ( x(2) > 4 || x(2) <= 6)
	state_cost = 4.4;
end

if (x(1) > 6) && ( x(2) > 6 )
	state_cost = 3;
end

