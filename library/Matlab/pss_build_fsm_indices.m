function xend_ind=pss_build_fsm_indices(params_symb,xend,xoru)

%
% xend_ind=pss_build_fsm_indices(params_symb,xend, xoru)
%
% Computes the coordinate-wise labels for a state given as a coordinate vector.
%   
% INPUTS: params_symb - data structure containing discretization information (see online documentation);
%         xend        - end-point(s) of a simulation (column vector(s));
%         xoru        - flag taking the value 'x' or 'u' to select labeling for x and u, respectively;
%
% OUTPUT: xend_ind    - labels for each coordinate of xend (column vector(s)).
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab May 2010.

neps=ceil(params_symb.epsilon/params_symb.eta);

if xoru=='x'
	xend_ind=round(xend/params_symb.eta)-params_symb.min(params_symb.xoind);
    if params_symb.deter==2 && (sum(xend_ind>params_symb.num(params_symb.xoind))>0 || sum(xend_ind<zeros(params_symb.n,1))>0)
        msgbox('State out of range.','Error: State OoR','error','replace');
%        error('State out of range.');
         xend_ind=NaN*ones(params_symb.n,1);
    elseif params_symb.deter==1 && (sum(xend_ind-params_symb.num(params_symb.xoind)>neps)>0 || sum(xend_ind<-neps*ones(params_symb.n,1))>0)
        msgbox('State out of range.','Error: State OoR','error','replace');
%        error('State out of range.');
         xend_ind=NaN*ones(params_symb.n,1);
    end
elseif xoru=='u'
	xend_ind=round(xend/params_symb.mu)-params_symb.min(params_symb.uind);
    if sum(xend_ind>params_symb.num(params_symb.uind))>0 || sum(xend_ind<zeros(params_symb.m,1))>0
        msgbox('Input out of range.','Error: Input OoR','error','replace');
 %       error('Input out of range.');
        xend_ind=NaN*ones(params_symb.m,1);
    end
end
