tau=200*1e-9;
lambda=[0.2,0.5,1,2,3,5,10,20,30,40,50,80]*1e6;
delta=0;
for j=1:8
    i=1;
    T2=[0];
    while sum(interval)<=0.1
        interval(i)=exprnd(1/lambda(j));
        delta=delta+interval(i);
        if delta>=tau
            T2(length(T2)+1)=sum(interval);
            delta=0;
        end
        i=i+1;
    end
    interval=[];
    delta=0;
    Tmatrix(2,j)={T2};
end
for j=9:12
    T2=[0];
    while max(T2)<=0.1
        T2(length(T2)+1)=max(T2)+tau+exprnd(1/lambda(j));
    end
    Tmatrix(2,j)={T2};
end