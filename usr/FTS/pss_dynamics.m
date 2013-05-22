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

%% INITIALIZATIONS

if xl==0 && ul ==0
    xend=1;
elseif xl==0 && ul==1
    xend=[0 2];
elseif xl==1 && ul==0
    xend=[1 2];
elseif xl==1 && ul==1
    xend=[3 1];        
elseif xl==2 && ul==0
    xend=1;            
elseif xl==2 && ul==1
    xend=[3 4];                
elseif xl==3 && ul==0
    xend=4 ;                    
elseif xl==3 && ul==1
    xend=[2 1];                        
elseif xl==4 && ul==0
    xend=4;                            
elseif xl==4 && ul==1
    xend=4;
end