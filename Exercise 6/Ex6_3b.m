%Ex 6 3b


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

time = 1:N

time_blocks = 6;
b_legnth = N/time_blocks;

I = sparse(eye(N));
Q_G = sparse(kron(I,Q));
R_G = sparse(kron(b_legnth*eye(time_blocks),R));
G = blkdiag(Q_G,R_G);

A_eq_I = sparse(eye(N*nx));
A_eq_A = sparse(kron(diag(ones(N-1,1),-1),-A));
one_block = sparse(kron(eye(time_blocks),ones(b_legnth,1)));    
Aeq_B = kron(one_block, -B);
A_eq = [A_eq_I + A_eq_A, Aeq_B];
b_eq = sparse([A*x0;zeros((N-1)*nx,1)]);

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
u_blocks= z(N*nx+1:nx*N+nu*time_blocks);
u = one_block*u_blocks;       

% Plot optimal trajectory
figure(5);
subplot(2,1,1);
plot([0,time],y,'-ko'); % Plot on 0 to N
grid('on');
ylabel('y_t');
subplot(2,1,2);
hold('on');
stairs(time-1,u,'k'); % Plot on 0 to N-1
plot(time-1,u,'ko');  % Plot on 0 to N-1
hold('off');
box('on');
ylim([-2,2]);
grid('on');
xlabel('t');
ylabel('u_t');