function iterate_fti(system,state,trace,m,n)
    trans=pss_fti(system,state(:,end));
    len=size(trans,2);
    if isempty(trans)
        trace
    else
        for i=1:len
           iterate_fti(system,trans(m+1:n+m,i),[trace trans(1:m,i)],m,n);
        end
    end