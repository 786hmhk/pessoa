function add_state = pss_target_set(x)

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


add_state=1;

if x(1)>=1.8 && x(1)<=2 && x(2)>=1 && x(2)<=4
    add_state=0;
end

if x(1)>=2.8 && x(1)<=3 && x(2)>=2 && x(2)<=5
    add_state=0;
end

if x(1)>=3.8 && x(1)<=4 && x(2)>=1 && x(2)<=4
    add_state=0;
end
