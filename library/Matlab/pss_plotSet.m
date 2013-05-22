function pss_plotSet(params_symb,p1)

%
% pss_plotSet(params_symb,p1)
%
% This function plots a point given its labels.
%   
% INPUTS: params_symb - data structure containing discretization information
%                       (see online documentation);
%         p1          - state labels (column vector).
%
% OUTPUT: none.
%
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab 2009.

pp1=params_symb.eta*(p1'+params_symb.min(params_symb.xoind)); 
plot(pp1(1),pp1(2),'*');
