% [X,Y] = meshgrid(-2:0.5:8);
% Z = -(3-0.4.*X).*X-(2-0.2.*Y).*Y;
% figure(1)
% contour(X, Y, Z,10)
x_d = 1.32;
kappa = 2.4;
g = 9.81;
h = 0.01;
t = 0:h:10;
y_0 = [2;0];
m = 200;

y_eul = zeros(2,length(t)); % 2 x length(t) matrix
y_b = zeros(2,length(t));
y_c = zeros(2,length(t));

y_eul(:,1) = y_0; %first column
y_b(:,1) = y_0;
y_c(:,1) = y_0;

f = @(y) [y(2); -g*(1-(x_d/y(1))^kappa)];

opt = optimset('Display','off','TolFun',1e-8); % Options for fsolve

for i = 1:(length(t)-1)
    y_eul(:,i+1) = y_eul(:,i) + h*feval(f,y_eul(:,i));
    
    r = @(yb_next) (y_b(:,i) + h*feval(f,yb_next)-yb_next);
    y_b(:,i+1) = fsolve(r,y_b(:,i),opt);
    
    rc = @(ycnext) (y_c(:,i) + h*feval(f,(ycnext + y_c(:,i))/2)-ycnext);
    y_c(:,i+1) = fsolve(rc, y_c(:,i),opt);
    
      
      
    
    
 

end

E_a = (m*g/(kappa-1))*x_d^(kappa).*y_eul(1,:).^(1-kappa)+m*g.*y_eul(1,:)+0.5*m.*y_eul(2,:).^2;
E_b = (m*g/(kappa-1))*x_d^(kappa).*y_b(1,:).^(kappa-1)+m*g.*y_b(1,:)+0.5*m.*y_b(2,:).^2;
E_c = (m*g/(kappa-1))*x_d^(kappa).*y_c(1,:).^(kappa-1)+m*g.*y_c(1,:)+0.5*m.*y_c(2,:).^2;
%plot(t,y_eul(1,:));
%plot(t,y_b(1,:));
%plot(t,y_c(1,:));

figure;
hold on
grid on
plot(t, E_a);

plot(t, E_b);

plot(t,E_c);
legend("explicit euler", "implicit", "implicit midpoint");

