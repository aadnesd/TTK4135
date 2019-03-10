%Ex 6 1d


A = [1  0.5;
     0  1  ];
b = [0.125  0.5]';


Q = diag([2 2]);
R = 2;


[~,P] = dlqr(A,b,Q/2,R/2);

I = [1 0; 0 1];
K = (1/(R/2))*b'*P*((I+b*(1/(R/2))*b'*P)\A);

eig_closed = eig(A-b*K); % both <1
