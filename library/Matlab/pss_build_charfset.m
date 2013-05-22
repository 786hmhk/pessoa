function [set_state_array nstates] = pss_build_charfset(params_symb,state_array,nbatch)

%
% [set_state_array nstates] = pss_build_charfset(params_symb,state_array,size_batch)
%
% Constructs a symbolic (ROBDD) representation for a set given a characteristic
% function defined in "pss_target_set.m" by the user.
%   
% INPUTS: params_symb   - data structure containing discretization information
%                         (see online documentation);
%         state_array   - array describing a hyperbox containing the target set
%						  described by the characteristic function;
%         size_batch    - size of computation batches to be used (see online
%                         documentation);
%
% OUTPUT: set_state_array - array containing states belonging to target set 
%			    (corresponding with the characteristic function);
%         nstates	  - total number of states in the target set.
%
% Anna Davitian <davitian@ee.ucla.edu>, UCLA CyPhyLab 2009.

% n is the number of states not in target set (initialized to entire state space)
n=0;
set_state_array = [];

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
	add_state = pss_target_set(user_state);

	if add_state == 1
		set_state_array = [set_state_array user_state];
	elseif add_state == 0
		n=n+1;
	else
		%beep;
		warndlg('Target set output not defined properly! Only use {0,1}.','Pessoa Warning')
		error('Target set output not defined properly! Only use {0,1}. Execution stopped.')
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



