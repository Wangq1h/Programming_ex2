sample1=[1.03,1.05,1.07];
sample2=[1.02,1.04,1.08];
tau=0.015;
N1=length(sample1);
N2=length(sample2);
B=sample1'*ones(1,N1);
C=ones(N2,1)*sample2;
delta=B-C;
[row,col]=find(abs(delta)<=tau);
m=find(abs(delta)<=tau);
row1=row(find(delta(m)<0));
col1=col(find(delta(m)>=0));
result=[sample1(row1),sample2(col1)]