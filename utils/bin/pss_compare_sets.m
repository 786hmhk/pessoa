function output = pss_compare_sets(set1_name, set2_name)

%
% output = pss_compare_sets(set1,set2)
%
% Tests if set2 is contained in set1 or if the intersection is non-empty.   
%   
% INPUTS: set1 - string containing the filename where the first set is stored;
%         set2 - string containing the filename where the second set is
%                stored (sets need to have been computed with the same eta);
%		
% OUTPUTS: output - flag with the result of the test:
%          output=0 when the sets have empty intersection;
%          output=1 when set2 is contained in set1;
%          output=2 when the sets have non-empty intersection;
%
% Anna Davitian <davitian@ee.ucla.edu>, UCLA CyPhyLab 2009.
% Manuel Mazo Jr. <mmazo@ee.ucla.edu>, UCLA CyPhyLab 2008.

disp('------------------------------------------------------------------------');
disp('    ');
disp('                             PESSOA Version 1.0                       '); 
disp('                  UCLA Cyber-Physical Systems Laboratory');
disp('                      http://www.cyphylab.ee.ucla.edu ');
disp('    ');
disp('------------------- Pessoa: Comparing Sets Initiated ------------------- ');

load(strcat(set2_name, '_symb'),'params_symb');

params_symb2=params_symb;

if strcmp(set1_name(end-2:end),'dom')
    load(strcat(set1_name(1:end-4), '_symb'),'params_symb');
else
    load(strcat(set1_name, '_symb'),'params_symb');
end

n=params_symb.n;
m=params_symb.m;

if(params_symb2.n==n && sum(params_symb2.min==params_symb.min)==(n+m) && sum(params_symb2.max==params_symb.max)==(n+m) && params_symb2.eta==params_symb.eta)
	output = pessoa_compare_sets(set1_name, set2_name);
else
	error('The sets must have the same semantics.')
end



