function variance(N,M)

S=0;
P=zeros(N,1);
t=0.5;
r=0;

for n=1:N
for k=1:M
     r=unifrnd(0,1);
    if r>t
        S=S+1;
    else
        ;
    end
      t=r;
end
P(n)=S;
S=0;
t=0.5;
end

M=mean(P);
V=std(P)^2;

M
V
