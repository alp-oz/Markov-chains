

function Rt=deinceps(p,a,n,m) %R Finding the R_n=s_n*s_{n-1}*...  values

g=declinis(a,0.8*a,p,m,m);


%determining W after declinis, the issue is it is unstable so we need to
%control the tail
% n is the number of terms (initial values) in the distribution, w_0, w_1
% etc
%num2str(ans, 10)
q=1-p;
t=0;
s=0;

A = zeros(n,1);
G = zeros(n,1);
D = zeros(n,1);

%H=0;

S = zeros(n,1);

X = zeros(n,1); %alpha
Y = zeros(n,1); %gamma
%M = zeros(N,2);

W = zeros(n,1);
Q = zeros(n,1);


R = zeros(n,1);
Rt= zeros(n,1);
Z = zeros(n,1);

%alpha_0 and gamma_0

X(1)=a; 
Y(1)=g;

% Determine zero Root of the poly A_0G_0D_0-w_0=0
f= @(w) (w+a*(1-w))*(w+q/(1+q)*(1+q*g)*(1-w))*(w+p*(1-w))-w;

syms y
c=sym2poly(f(y));
v=roots(c);

v=v(v>=0);
v=sort(v);
%v;
w=v(1);

%Initial values 

W(1)=w; % w_0
R(1)=1-w; % dim 2, dim is the size of the matrices A is 2x2 etc.
S(1)=1-w; 

%functions A_i the row sum of A matrix
%X(1)=(a+S(1)-1)/(q*S(1));  %dim 3, alpha_1
%Y(1)=(g+S(1)-1)/(q^2*S(1)); %dim 3, gamma_1

A(1)=w+a*(1-w); %A_0, the first row of A
A(2)=1-A(1);

G(1)=w+q/(1+q)*(1+q*g)*(1-w);
G(2)=1-G(1);

D(1)=w+p*(1-w);
D(2)=1-D(1);

for i=2:n-1     % W(n)=w_(n-1) dim: n --> (n+1) the matrix n by n. At the first step (n=2), we already have w_0 and split the remainder to determine w(1)

    for j=1:i %i(i+1)/2 terms, which makes sense. For Z_1(Z_2) we need 3, for Z_2(Z(3)) we need 6
        for k=1:i+1-j
            Z(i) = Z(i) + A(j)*G(k)*D(i+2-j-k);
        end
    end

 

%S-formula, I hope there is no mistake here

S(i)=1/(1+G(1)*(p/q*D(1)-q^(2*i-1)*A(1)))*(1+D(1)*((1-X(i-1))/q*G(1)+q^(i-1)*(1-Y(i-1))/(1+q)*A(1))-Z(i)/R(i-1));


% Update alpha and gamma
X(i)=(X(i-1)+S(i)-1)/(q*S(i));
Y(i)=(Y(i-1)+S(i)-1)/(q^2*S(i));

%gamma to compensate for sensitivity, imperfect part of the code. If it
%stabilizes, stop

if i>=3 && abs(Y(i)-Y(i-1)) >=abs(Y(i-1)-Y(i-2))
    Y(i)=(Y(i-1)+Y(i-2))/2;
    t=1;
else
end


if t==1
if i>=3 && abs(X(i)-X(i-1)) >= abs(X(i-1)-X(i-2))
    X(i)=(X(i-1)+X(i-2))/2;
    s=1;
else
end
else 
end

if s==1
if i>=3 && abs(S(i)-S(i-1)) >=abs(S(i-1)-S(i-2))
    S(i)=(S(i-1)+S(i-2))/2;
else
end
end

W(i)=(1-S(i))*R(i-1);

R(i)= R(i-1)-W(i);



%A(i+1)= (1-X(i))/(1-q*X(i))*A(i);
A(i+1)= (1-X(i))*R(i);
A(i) = A(i)-A(i+1);

G(i+1) = R(i)*q^(i-1)*(1-q^2*Y(i))/(1+q);
G(i) = G(i) - G(i+1);

D(i+1) = R(i)*q^(2*i-1);
D(i) = D(i)-D(i+1);

  % for j=1:i %i(i+1)/2 terms, which makes sense. For Z_1(Z_2) we need 3, for Z_2(Z(3)) we need 6
  %      for k=1:i+1-j
            Q(i) = Q(i) + A(j)*G(k)*D(i+2-j-k);
  %      end
  %  end
 %Q(i)=Q(i)-W(i);
   

end

% False, find a new formula...H=(S(1)*(1-q*a)-p*S(1)*S(2))/p

%False... K=-q*(1-Y(1))*S(1)/(p*(1+q));

%for i=1:n
%    K=K+R(i);
%end

for m=1:n
    if R(m)>=0 
    Rt(m)=R(m);
    else
    Rt(m) = 0;
    end
end


for m=1:n
    if R(m)>=0 
    Rt(m)=R(m);
    else
    Rt(m) = 0;
    end
end

%disp(Q);
