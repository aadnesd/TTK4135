x = linspace(-1,6);
y = linspace(-1,6);
[X,Y] = meshgrid(x,y);
Z = -3*X -2*Y;

figure 
contour(X,Y,Z,10,'linewidth',1);


