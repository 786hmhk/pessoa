function pss_design(sys_name,target_name,contr_name,type,N_sims,verbose)

%	pss_design(sys_name,target_name,contr_name,type,N_sims,verbose)
%
%	Synthesizes a controller for the abstraction constructed, solving one of the following
%	problems: Reachability, Safety, Reach and Stay, Reach while Stay. The
%	controller is saved in the form of a BDD.
%   
%       INPUTS: sys_name	- name of the system abstraction to be controlled. 
%               target_name	- name of the target set(s) to be used: target_name{1}=set_name. 
%                           For Reach&Stay while Safe: target_name{1}=Reach_Set_name; 
%                                                      target_name{2}=Safe_Set_name.
%               contr_name	- name for the controller to be synthesized.
%               type		- 1: Safety, 2: Reachability, 3:Reach and Stay, 4:Reach while Stay.
%               N_sims		- size of batches to be used (See Online Documentation for more info).
%               verbose     - verbose level: 0 (minimal), 1 (extra info), 2 (1+plots, just 2D systems), 3 (debug mode).
%
%       OUTPUT: none
%
%   Manuel Mazo Jr. <mmazo@ee.ucla.edu>, CyPhyLab-UCLA 2009.

load(strcat(sys_name,'_symb'),'params_symb');

domain_name=strcat(contr_name,'_dom');

nbatch=N_sims;
totloops=ceil(prod(params_symb.num+ones(size(params_symb.num)))/nbatch);

d=size(target_name);

if type<=3 && d(2)~=1
    error('Just one target_name needed.');
elseif type==4 && d(2)~=2
    error('Two target_names needed.');
end

% type_flag is 0 normally; 1 for constraint set; 2 for target set - only used for printing purposes for safety only
type_flag = 0;

if type<3
    pessoa_design(sys_name,target_name{1},contr_name,domain_name,type,type_flag,verbose);
elseif type==3
    pessoa_design(sys_name,target_name{1},'TMPcontS','TMPsafeDom',1,type_flag,verbose);
    pessoa_design(sys_name,'TMPsafeDom','TMPcontR','TMPreachDom',2,type_flag,verbose);
    pessoa_combine('TMPcontS','TMPcontR',contr_name,domain_name,0,type_flag,verbose);
elseif type==4
    %First solve for the Safety in the Safe set
    pessoa_design(sys_name,target_name{2},'TMPcontSafe','TMPsafeDom',1,1,verbose);
    %Then use this controller to design a Reach and Stay controller
    pessoa_design('TMPcontSafe',target_name{1},'TMPcontS','TMPsafeDom',1,2,verbose);
    pessoa_design('TMPcontSafe','TMPsafeDom','TMPcontR','TMPreachDom',2,type_flag,verbose);
    pessoa_combine('TMPcontS','TMPcontR',contr_name,domain_name,0,type_flag,verbose);
elseif type==5 || type==6
    pessoa_design(sys_name,target_name{1},contr_name,domain_name,type,type_flag,verbose);
else
    error('Error: type must be 1 (Safety), 2 (Reachability) or 3 (Reach and Stay)');
end

filename=strcat(contr_name,'_symb');
save(filename, 'params_symb');

if type>=3
    !rm TMP*.bdd
end;
