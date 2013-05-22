function u=pss_controller(contr_name,uo,x)

%
% u=pss_controller(controller,uo,x)
%
% Lazy determinization of a symbolic controller. Computes the input with
% the shortest ROBDD description that can be applied at the state x. Reuses
% the previous input if it can be applied at x.
%
% INPUTS: controller - name of the file where the controller is stored;
%         uo         - previous input;
%         x          - state.
%
% OUTPUT: u - input that can be applied at state x.
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab 2009.

load(strcat(contr_name,'_symb'), 'params_symb');

xl=pss_build_fsm_indices(params_symb,x,'x');
ul=pss_build_fsm_indices(params_symb,uo,'u');
ul=pessoa_controller(contr_name,ul,xl);
ul=ul';
if isnan(ul)
    u=ul;
else
    u=params_symb.mu*(ul+params_symb.min(params_symb.uind));
end
