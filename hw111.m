clear
result=[];
tau=1e-5;
load("修正后的第一时序.mat")
load("时序2.mat")
N=9000;
%%
for i=1:N
    Tadj1=T1adj(find(T1adj<=i*36000/N & T1adj>(i-1)*36000/N));
    Tadj2=T2(find(T2<=i*36000/N & T2>(i-1)*36000/N));
    N1=length(Tadj1);
    N2=length(Tadj2);
    B=Tadj1'*ones(1,N2);
    C=ones(N1,1)*Tadj2;
     delta=B-C;
    [row,col]=find(abs(delta)<=tau);
     m=find(abs(delta)<=tau);
    row1=row(find(delta(m)<0));
    col1=col(find(delta(m)>=0));
    result=[result,Tadj1(row1),Tadj2(col1)];
end
%%
result=sort(result);
A=diag(-ones(1,length(result)-1),-1)+eye(length(result));
delta2=A*result';
histogram(delta2)
legend
title("二重偶然符合时间差的分布")
xlabel("时间(s)")
ylabel("频次")