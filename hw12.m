lambda=1/3;
p1=0;
p2=0;
T1=[];
T2=[]; 
W=[];
N=10000000;
for i=1:N
    interval1=exprnd(1/lambda);
    interval2=interval1;
    cho=rand(1);
    W(i)=interval1;
    while cho<=0.75
        interval2=interval2+exprnd(1/lambda);
        cho=rand(1);
    end
    if interval2>interval1
        cost=26;
    else
        cost=16;
    end
    T1(i)=interval1+cost;
    T2(i)=interval2+16;
    if abs(T1(i)-T2(i))==0
        p1=p1+1;
    end
    if T2(i)-T1(i)>0
        p2=p2+1;
    end
end
aveT1=mean(T1);
aveT2=mean(T2);
p1=p1/N;
p2=p2/N;
aveW=mean(W);
save("101.mat","T1")
save("102.mat","T2")
save("103.mat","p1")
save("104.mat","p2")
save("105.mat","aveT1")
save("106.mat","aveT2")