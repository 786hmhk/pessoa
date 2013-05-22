function plot_domain(contr_name)

load(strcat(contr_name,'_symb'), 'params_symb');

N=(params_symb.num(1)+1)*(params_symb.num(1)+2);

    figure(1);
    xl=[0;0];
    for j=1:N
        ul=pessoa_control_nowarn(contr_name,[0;0],xl);
		clear mex;
        ul0=pessoa_control_nowarn(contr_name,[0;0],xl);
        ul1=pessoa_control_nowarn(contr_name,[0;1],xl);        
		clear mex;
        if(sum(isnan(ul))==0 && ul0(2)==0 && ul1(2)==0)
           x1=params_symb.eta*(xl(1)+params_symb.min(1));
           x2=params_symb.eta*(xl(2)+params_symb.min(2));        
           figure(1);hold on; plot(x1,x2,'r.');
        elseif(sum(isnan(ul))==0 && ul0(2)==1 && ul1(2)==1)
           x1=params_symb.eta*(xl(1)+params_symb.min(1));
           x2=params_symb.eta*(xl(2)+params_symb.min(2));        
           figure(1);hold on; plot(x1,x2,'b.'); 
        elseif(sum(isnan(ul))==0 && ul0(2)==0 && ul1(2)==1)
           x1=params_symb.eta*(xl(1)+params_symb.min(1));
           x2=params_symb.eta*(xl(2)+params_symb.min(2));        
           figure(1);hold on; plot(x1,x2,'m.');            
        end
        
        k=1;
        while (k<=2 && xl(k)>=params_symb.num(k))
            xl(k)=0;
            k=k+1;
        end
        if k>2
                break;
        else
            xl(k)=xl(k)+1;
        end    
    end
    