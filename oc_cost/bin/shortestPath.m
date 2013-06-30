function shortestPath()

clear all

%% Testbench files: Delete them. :)
load('tb_sys_info.txt')
params_symb.n       = 1;
params_symb.m       = 1;
params_symb.nume(1) = tb_sys_info(1);
params_symb.nume(2) = tb_sys_info(2);
save('test_symb.mat', 'params_symb');
%% 

load('test_symb.mat');

shortestPathMex();



% load(strcat(sys_name,'_symb'),'params_symb');
% 
% domain_name=strcat(contr_name,'_dom');
% 
% nbatch=N_sims;
% totloops=ceil(prod(params_symb.num+ones(size(params_symb.num)))/nbatch);
% 
% d=size(target_name);
% 
% if type<=3 && d(2)~=1
%     error('Just one target_name needed.');
% elseif type==4 && d(2)~=2
%     error('Two target_names needed.');
% end
% 
% % type_flag is 0 normally; 1 for constraint set; 2 for target set - only used for printing purposes for safety only
% type_flag = 0;
% 
% if type<3
%     pessoa_design(sys_name,target_name{1},contr_name,domain_name,type,type_flag,verbose);
% elseif type==3
%     pessoa_design(sys_name,target_name{1},'TMPcontS','TMPsafeDom',1,type_flag,verbose);
%     pessoa_design(sys_name,'TMPsafeDom','TMPcontR','TMPreachDom',2,type_flag,verbose);
%     pessoa_combine('TMPcontS','TMPcontR',contr_name,domain_name,0,type_flag,verbose);
% elseif type==4
%     %First solve for the Safety in the Safe set
%     pessoa_design(sys_name,target_name{2},'TMPcontSafe','TMPsafeDom',1,1,verbose);
%     %Then use this controller to design a Reach and Stay controller
%     pessoa_design('TMPcontSafe',target_name{1},'TMPcontS','TMPsafeDom',1,2,verbose);
%     pessoa_design('TMPcontSafe','TMPsafeDom','TMPcontR','TMPreachDom',2,type_flag,verbose);
%     pessoa_combine('TMPcontS','TMPcontR',contr_name,domain_name,0,type_flag,verbose);
% elseif type==5 || type==6
%     pessoa_design(sys_name,target_name{1},contr_name,domain_name,type,type_flag,verbose);
% else
%     error('Error: type must be 1 (Safety), 2 (Reachability) or 3 (Reach and Stay)');
% end
% 
% filename=strcat(contr_name,'_symb');
% save(filename, 'params_symb');
% 
% if type>=3
%     !rm TMP*.bdd
% end;
