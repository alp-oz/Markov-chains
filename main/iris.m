

function T=iris(p,N,f)

%N is the number of parts of 0,1; a=1/n to 1 - how finer the curve is
%partitioned
%n is the first n values of w - how fine each estimate is

q=1-p;

R=zeros(N*f,N*f+10);
T=zeros(N-1,2);

u=zeros(N*f,1);
b=zeros(N*f,1);

U=zeros(N+1,1);
B=zeros(N+1,1);

t=0;
s=0;

for i=1:N*f
    R(i,1:10+ceil(i/5))=deinceps(p,(1-i/(N*f))/(1-q*(i/(N*f))),10+ceil(i/5),30);
end

for i=1:N*f
   % u(i)=R(i,1)-R(i,2)+q*R(i,3);
   u(i)=R(i,1)-R(i,2);
   b(i)=R(i,1);
   n=10+i;
    for j=2:n
     u(i)=u(i)+q^(j-1)*R(i,j);
     b(i)=b(i)+(1-q^(2*j-3))*R(i,j);
    end
     if R(i,n)~=0
    s=R(i,n)/R(i,n-1);
    u(i)=u(i)+(q^(n-1)*R(i,n)*(1/(q*s)));
    b(i)=b(i)+R(i,n)*(1/(1-s))-(q^(2*n-3)*R(i,n)*(1/(q^2*s)));
    else
     end
if u(i)>=s
    t=t+1;
    s=t*(q/p)*(1/N);
    U(t)=u(i);
    B(t)=b(i);
else
end
end

LU=U(1:end-2);
LB=B(1:end-2);

T(:,1)=LU;
T(:,2)=LB;
