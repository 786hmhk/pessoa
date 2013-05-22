function next_label = pss_increase_label(params_symb, current_label, flag)

%
% next_label = pss_increase_label(params_symb, current_label, flag)
%
% Increases current_label. Based on the flag, it either increases from left to
% right (flag=0), or bottom to top (flag=1).
%   
% INPUTS: params_symb   - data structure containing discretization information
%                         (see online documentation);
%         current_label - label to be increased, given in vector form [xini;u];
%         flag          - determines the direction in which the labels are to
%                         be increased:
%                         flag=0 increases labels from left to right;
%                         flag=1 increases labels from bottom to top.
%
% OUTPUT: next_label    - returns the increased label in vector form [x; y;z; ...].
%
% Anna Davitian <davitian@ee.ucla.edu> UCLA CyPhyLab 2009.
% Manuel Mazo <mmazo@ee.ucla.edu> UCLA CyPhyLab 2009.


if nargin <= 2
	flag = 0;
end

for i=1:params_symb.n
	label_min(i) = params_symb.min(i);
	label_max(i) = params_symb.max(i);
end

x = current_label;

% This is the default case (flag is 0), increasing left to right
if flag == 0
	for j=1:params_symb.n
		if x(j)<label_max(j)
			x(j)=x(j)+1;
			break;
		else
			x(j)=label_min(j);
		end
	end
end

% This is increasing upwards when flag is 1
if flag == 1
	for j=params_symb.n:-1:1
		if x(j)<label_max(j)
			x(j)=x(j)+1;
			break;
		else
			x(j)=label_min(j);

		end
	end
end


next_label = x;

