clear
%tau=0.000001;
lambda1=100;
lambda2=200;
T1=[exprnd(1/lambda1)];
T2=[exprnd(1/lambda2)];
i=2;
while min(max(T1),max(T2))<=10*3600
    T1(i)=T1(i-1)+exprnd(1/lambda1);
    T2(i)=T2(i-1)+exprnd(1/lambda2);
    i=i+1;
end
