function add_state = pss_target_set(x)

%
%   add_state = pss_target_set(x)
%
%   This function allows the user to define a characteristic function
%   to define a target set. The function must return 1 whenever the input
%   state belongs to the target set and 0 otherwise.
%   
%       INPUTS: x - state under consideration [x1;x2;...]
%
%       OUTPUT: add_sate - 1: state belongs to the target set
%						   0: states does not belong to the target set
%
%   Anna Davitian <davitian@ee.ucla.edu>,   CyPhyLab-UCLA 2009


% The user is free to define any characteristic function

if x(1)>=-0.06 && x(1)<=0.06 && x(2)>=-0.3 && x(2)<=0.3
 	add_state=1;
else
 	add_state=0;
end;
