%Exercise 5_d

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

%Getting the variables
z = KKT_sol(1:N*(nx+nu));
y = [x0(3);z(nx:nx:N*nx)];
u = z(N*nx+1:nx*N+nu*N);

time = 1:N;

figure(1)
subplot(2,1,1);
plot([0,time],y,'-ko');
legend('y');
subplot(2,1,2);
plot(time-1,u,'-ko');
legend('u');

x = quadprog(G,[],[],[],A_eq,b_eq);
figure(2);
plot(x,'-ko');

