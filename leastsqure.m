
% Solving least square problem using optimization algorithm
% Author : Soundarya S
% Uses bracketing line search line search to determine step lenth alpha_k
% We use following steps to write the code

% Step-1: Define r_j(x) symbolically
% Step-2: Define the objective function f(x)
% Step-3: Define the Jacobian J(x) symobolically
% Step-4: Initialize and use the while loop to run the iteration
% Step-5: Obtain the optimizer x^*
% Step-6: Determine the best fit function phi(t)= t^2x_1^*+tx_2^*+tx_3*
% Step-7: Plot the function phi and data points (t_i,y_) in the same plotting

clear all;
close all;
format long

N=3; % number of variables (edit this according to objective function)
xi=sym(zeros(1,N));
syms t;
for i=1:N
    syms("x"+i);
    xi(1,i) = ("x"+i);
end
tval=[1 2 3 4 5];
yval=[2.1 3.9 6.2 8.1 10.1];
M=length(yval);
f=0;
phi= x1*(t^2)+ x2*t+x3;
for i=1:M
    r(i)=yval(i)-subs(phi,t,tval(i));
    f=f+0.5*r(i)^2;
end

% step size 
alpha=1;
c = 10^-4;
rho = 0.3;
x = [5 5 5]; % Initial Guess

I = eye(N);
e = 10^(-5); % Convergence Criteria
k = 1; % Iteration Counter

% Gradient Computation:
for i=1:N 
    df_dx(i)= diff(f, "x"+i);
end
G = subs(df_dx, xi, x);

for i=1:M
    for j=1:N
      J(i,j)=diff(r(i),"x"+j); 
    end
end

Jk=subs(J,xi,x);
rk=subs(r,xi,x);
Pk=-(Jk'*Jk)\(Jk'*rk');

% display table
fprintf('k: %d\t', k);
for i=1:N
   fprintf('x%d: %d\t',i, x(i)); 
end
fprintf('f(xk): %d\t\t', subs(f,[xi],[x]));
fprintf('||∇f(xk)||: %d\t', norm(G)); 
fprintf('\n');

while norm(G)>=e
    % compute step length alpha_k
    alpha=1;
    lh = subs(f,xi, x + alpha * Pk');
    rh = subs(f,xi,x)+ c * alpha*G'*Pk';

    while lh>rh
        alpha = rho*alpha;
        lh = subs(f,xi,x+alpha*Pk');
        rh = subs(f,xi,x)+c*alpha*G'*Pk';
    end

    for i=1:N 
      x(i)= x(i) + alpha * Pk(i); % Updated xk value
    end

    Jk=subs(J,xi,x);
    rk=subs(r,xi,x);
    Pk=-(Jk'*Jk)\(Jk'*rk');
    G = subs(df_dx, xi, x);
    k=k+1;


    % display table
    fprintf('k: %d\t', k);
    for i=1:N
       fprintf('x%d: %d\t\t\t',i, x(i)); 
    end
    fprintf('f(xk): %d\t\t', subs(f,[xi],[x]));
    fprintf('||∇f(xk)||: %d\t', norm(G)); 
    fprintf('\n');
    
end
%plotting
phi_t=subs(phi,xi,x);
scatter(tval,yval);
hold on
fplot(phi_t,[0 7]);
hold off
xlabel('time');