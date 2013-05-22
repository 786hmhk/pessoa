function xend=pss_dynamics(params_symb, params_cont, xl, ul)

%
%   xend = pss_dynamics(params_symb, params_cont, xl, ul)
%
%   This function computes the transitions in labels (as opposed to coordinates) directly. 
%	In general after the simulation step one would need to call "pss_build_fsm_indices.m"
%   at the end to return labels.
%   
%       INPUTS: xl - state labels.
%           	ul - input labels.
%
%       OUTPUT: xend	- ending state(s) in labels.
%
%   Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA 2009.

% Transform labels into states (I'll probably make a function for this like the one at the end)
u=params_symb.mu*(ul+params_symb.min(params_symb.uind));
x=params_symb.eta*(xl+params_symb.min(params_symb.xoind));


% Here computations on states
xend=[];


% Transform the vector xend into labels (I'll make it soon so that
% pss_build_fsm_indices operates directly on vectors so that this loop is
% inside it

for i=1:size(xend,2)
    xend(:,i)=pss_build_fsm_indices(params_symb,xend(:,i),'x');
end
