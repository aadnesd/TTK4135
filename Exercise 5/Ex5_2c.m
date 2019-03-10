%Ex 5 2c

A = [0     0      0    ;
     0     0      1    ;
     0.1  -0.79   1.78];
B = [1 0 0.1]';
C = [0 0 1];

x0 = [0;0;1];

Q = [0 0 0; 0 0 0; 0 0 2];
R = 2;

N= 30;
nx = 3;
nu = 1;

I = sparse(eye(N));
Q_G = sparse(kron(I,Q));
R_G = sparse(kron(I,R));
G = blkdiag(Q_G,R_G);

A_eq_I = sparse(eye(N*nx));
A_eq_A = sparse(kron(diag(ones(N-1,1),-1),-A));
A_eq_B = sparse(kron(I,-B));

A_eq = [A_eq_I + A_eq_A, A_eq_B];

%b_eq = sparse([A*x0;zeros((N-1)*nx,1)]); Changes each time



time = 1:N;
%Constraint, x not constrained so "set" to inf
x_con_l = -Inf(nx*N,1);
x_con_u = Inf(nx*N,1);
u_con_l = -ones(nu*N,1);
u_con_u = ones(nu*N,1);

lb = [x_con_l;u_con_l];
ub = [x_con_u;u_con_u];

A_real = [0 0 0; 0 0 1; 0.1 -0.855 1.85];
B_real = [1;0;0];
C_real = [0 0 1];


%Mpc
%Need to get new u and x each time
%Initialize
u = zeros(nu,N);
x = zeros(nx,N+1);
x(:,1) = x0;

b_eq = zeros(N*nx,1);

for t = 1:N
    b_eq(1:nx) = A*x(:,t);
    
    [z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],A_eq,b_eq,lb,ub,[]);
    
    u_new = z(N*nx+1:N*nx+N*nu);
    u(t) = u_new(1);
    
    x(:,t+1 ) = A_real*x(:,t)+B_real*u(t);
end;

y = C_real*x;

figure(1);
subplot(2,1,1);
plot([0,time],y,'-ko');
ylim([-0.5, 5]);
ylabel('y_t');
subplot(2,1,2);
plot(time-1,u,'-ko');
ylim([-1.5, 0.5]);
xlabel('t');
ylabel('u_t');