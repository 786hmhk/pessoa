function [set_state_array state_cost_array nstates] = pss_build_charfcost(params_symb,state_array,nbatch)

%
%  function [set_state_array state_cost_array nstates] = pss_build_charfcost(params_symb,state_array,nbatch)
%
%  Constructs a symbolic (ADD) representation for a states in the state space of the system
%  given a characteristic function defined in "pss_cost_def.m" by the user.
%   
%      INPUTS: params_symb   - data structure containing discretization information
%                              (see online documentation);
%              state_array   - array describing a hyperbox containing the whole state space.
%              size_batch    - size of computation batches to be used (see online
%                              documentation);
%
%      OUTPUT: set_state_array  - array containing states belonging to target set 
%	                          (corresponding with the characteristic function);
% 	       state_cost_array - array containing the costs of each state of the state space.
%              nstates          - total number of states in the state space.
%
%  Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>, DCSC - TU Delft, 2013

% n is the number of states not in target set (initialized to entire state space)
n=0;
set_state_array   = [];
state_cost_array  = [];

% Convert back to coordinates on the grid (user has no knowledge of labels)
for m=1:nbatch
	for s=1:params_symb.n
		array1(s,1) = state_array(s,m);
	end
	xarray=params_symb.eta*(array1+params_symb.min(params_symb.xoind));
	for t=1:params_symb.n	
		state(t,m)=xarray(t,1);
	end
end

% Call characteristic function to see if state belongs to the target set
for i=1:nbatch
	for r=1:params_symb.n
		user_state(r,1) = state(r,i);
	end
	
	state_cost = pss_cost_def(user_state);

	if state_cost ~= Inf
		set_state_array  = [set_state_array user_state];
		state_cost_array = [state_cost_array state_cost];
	else
		n=n+1;
	end
end

nstates = nbatch - n;

% Convert back to FSM indices
for k=1:nstates
	for l=1:params_symb.n
		array2(l,1) = set_state_array(l,k);
	end
	xlarray=round(array2/params_symb.eta)-params_symb.min(params_symb.xoind);
	for j=1:params_symb.n	
		set_state_array(j,k)=xlarray(j,1);
	end
end

