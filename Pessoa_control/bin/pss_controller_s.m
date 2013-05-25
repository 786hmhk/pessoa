function [sys,ul0,str,ts]=pss_controller_s(t,ul,x,flag,params_symb,contr_name,uini) 

%   PSS_CONTROLLER_S
%   Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA, 2009.

switch flag,
 
  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,                                                
    [sys,ul0,str,ts,simStateCompliance] = mdlInitializeSizes(params_symb,uini);    
    
  %%%%%%%%%%  
  % Update %
  %%%%%%%%%%
  case 2,                                               
    sys = mdlUpdate(t,ul,x,contr_name, params_symb);
    
  %%%%%%%%%%
  % Output %
  %%%%%%%%%%
  case 3,                                               
    sys = mdlOutputs(t,ul,x,contr_name, params_symb);    

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,                                               
    sys = mdlTerminate();

  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

%end sfundsc2

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,ul0,str,ts,simStateCompliance]=mdlInitializeSizes(params_symb,uini)
clear mex;
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = params_symb.m;
sizes.NumOutputs     = params_symb.m+1;%2*params_symb.m+params_symb.n+1;
sizes.NumInputs      = params_symb.n;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

ul0  = pss_build_fsm_indices(params_symb,uini,'u');

str = [];
ts  = [params_symb.tau 0]; % Sample period of "tau" seconds

% speicfy that the simState for this s-function is same as the default
simStateCompliance = 'DefaultSimState';

% end mdlInitializeSizes

%
%=======================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=======================================================================
%
function sys = mdlUpdate(t,ul,x,contr_name,params_symb)
% Update the internal state of the controller
%clear mex;
xl=pss_build_fsm_indices(params_symb,x,'x');
sys=pessoa_controller(contr_name,ul,xl);
    
%end mdlUpdate

%
%=======================================================================
% mdlOutputs
% Return the output vector for the S-function
%=======================================================================
%
function sys = mdlOutputs(t,ul,x,contr_name,params_symb)
stop=0;

xl=pss_build_fsm_indices(params_symb,x,'x');
if sum(isnan(xl))>0
    xl=-ones(params_symb.n,1);
end

ul=pessoa_control_nowarn(contr_name,ul,xl);
ul=ul';

uout=params_symb.mu*(ul+params_symb.min(params_symb.uind));

if (isnan(uout))
    ul=zeros(params_symb.m,1);
    uout=ul;
    stop=1;
end;

sys=[uout;stop];
%sys=[uout;xl;ul;stop]; %input, labels-states, labels-inputs, stop
%end mdlOutputs

function sys = mdlTerminate()
clear mex;
sys=[];