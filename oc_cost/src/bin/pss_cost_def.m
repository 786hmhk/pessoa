function state_cost = pss_cost_def(x)

%
%  function state_cost = pss_cost_def(x)
%
%  This function allows the user to define any characteristic function
%  in any fashion to define the cost of each state of the system's state
%  space, as long as the state x belongs to it. Returns the cost of a
%  valid state, inf otherwise.
%   
%      INPUTS: x          - state to be considered [x1;x2;...]
%
%      OUTPUT: state_cost - Returns the cost of the state. If the state 
%			    passed belongs to characteristic function it
%			    it ruturns its value, inf otherwise
%
%  Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>, DCSC - TU Delft, 2013

% The user is free to define any characteristic function


state_cost = Inf;


 
if (x(1) >= 0 && x(1) <= 2) && ( x(2) >= 0 && x(2) <= 2)
	state_cost = 8;
else
    state_cost = 1;
end
 
if (x(1) > 2 && x(1) <= 4) && ( x(2) > 2 && x(2) <= 4)
	state_cost = 6;
end
 
if (x(1) > 4 && x(1) <= 6) && ( x(2) > 4 && x(2) <= 6)
	state_cost = 4.4;
end

if (x(1) > 6) && ( x(2) > 6 )
	state_cost = 3;
end
 

