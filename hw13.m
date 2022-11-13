clear
tau=200*1e-9;
lambda=[0.2,0.5,1,2,3,5,10,20,30,40,50,80]*1e6;
interval=[];
delta=0;
Tmatrix=cell(2,length(lambda));
for j=1:12
    i=1;
    T1=[0];
    T2=[0];
    while sum(interval)<=0.1
        interval(i)=exprnd(1/lambda(j));
         if interval(i)>=tau
             T1(length(T1)+1)=sum(interval);
         end
      delta=delta+interval(i);
         if delta>=tau
             T2(length(T2)+1)=sum(interval);
              delta=0;
         else
             continue;
         end
      i=i+1;
    end
    interval=[];
    delta=0;
    Tmatrix(1,j)={T1};
    Tmatrix(2,j)={T2};
end
save("两种系统的记录","Tmatrix")
%%
close all
lambda=[0.2,0.5,1,2,3,5,10,20,30,40,50,80]*1e6;
lambdaa=[];
for i=1:12
   lambdaa(i)=length(Tmatrix{1,i})/max(Tmatrix{1,i});
   lambdab(i)=length(Tmatrix{2,i})/max(Tmatrix{2,i});
end
scatter(lambda,lambdaa)
hold on
scatter(lambda,lambdab)
hold on
lambda=1:80*1e6;
yy=lambda./(1+lambda.*tau);
plot(lambda,yy)
hold on
yyy=lambda./exp(lambda.*tau );
plot(lambda,yyy)
hold off
legend