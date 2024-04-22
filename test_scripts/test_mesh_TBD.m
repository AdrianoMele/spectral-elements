% Define reference element
xr = [1 0 0 1];
yr = [0 0 1 1];


% Define domain
NN = 20;
th = linspace(0,2*pi,NN*4)';
xd = cos(th); yd = sin(th);


% Find map between domains
g1.x = xd(1:NN); g1.y = yd(1:NN);
g2.x = xd(NN:2*NN); g2.y = yd(NN:2*NN);
g3.x = xd(2*NN:3*NN); g3.y = yd(2*NN:3*NN);
g4.x = xd(3*NN:4*NN); g4.y = yd(3*NN:4*NN);

% Plot domain
figure
hold on
plot(g1.x,g1.y, 'o-')
plot(g2.x,g2.y, 'o-')
plot(g3.x,g3.y, 'o-')
plot(g4.x,g4.y, 'o-')
axis equal