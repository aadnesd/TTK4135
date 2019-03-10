%Ex 5 f


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

b_eq = sparse([A*x0;zeros((N-1)*nx,1)]);

KKT_zero = sparse(zeros(N*nx));
KKT_zero_vec = sparse(zeros(N*(nx+nu),1));

KKT = [G -A_eq';
       A_eq KKT_zero];
KKT_vec = [KKT_zero_vec;b_eq];

KKT_sol = KKT\KKT_vec;

time = 1:N;
%Constraint, x not constrained so "set" to inf
x_con_l = -Inf(nx*N,1);
x_con_u = Inf(nx*N,1);
u_con_l = -ones(nu*N,1);
u_con_u = ones(nu*N,1);

lb = [x_con_l;u_con_l];
ub = [x_con_u;u_con_u];  

%Getting the variables

[z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],A_eq,b_eq,lb,ub,[]);

y = [x0(3);z(nx:nx:N*nx)];
u= z(N*nx+1:nx*N+nu*N);


figure(1)
subplot(2,1,1);
plot([0,time],y,'-ko');
ylim([-1,4]);
legend('y');
subplot(2,1,2);
plot(time-1,u,'-ko');
legend('u');