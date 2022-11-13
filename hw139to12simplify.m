tau=200*1e-9;
lambda=[0.2,0.5,1,2,3,5,10,20,30,40,50,80]*1e6;
interval=[];
for j=9:12
    i=1;
    T1=[0];
    T2=[0];
    while sum(interval)<=0.001
        interval(i)=exprnd(1/lambda(j));
        if interval(i)>tau
            T1(length(T1)+1)=sum(interval);
        end
        i=i+1;
    end
    TT1=[];
    for jj=1:100
        TT1=[TT1,T1+(jj-1)*0.001];
    end
    interval=[];
    Tmatrix(1,j)={TT1};
end
