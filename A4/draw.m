function draw(fnum,xmin,xmax,dx,ymin,ymax,dy)

drawMin = false;
drawAreas = false;

[x1,x2] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);
f = zeros(size(x1));
for i = 1 : size(x1,1)
	for j = 1 : size(x1,2)
		x = [x1(i,j),x2(i,j)]';
		f(i,j) = fun(x,fnum);
	end
end
clf;
colormap('default');
subplot(1,2,1);
surf(x1,x2,f);
hold on;
grid on;
axis equal;
axis([min(min(x1)),max(max(x1)),min(min(x2)),max(max(x2)),min(min(f)),max(max(f))]);
xlabel('x_1');
ylabel('x_2');
zlabel('f');

% Optionally draw the minimum point
if drawMin
	% Optimizer options
	optUnc = optimoptions(@fminunc,...
		'SpecifyObjectiveGradient',true,...
		'CheckGradients',false,...
		'Display','off');
	xinit = [-1,-1]';
	xmin = fminunc(@(x)fun(x,fnum),xinit,optUnc);
	plot3(xmin(1),xmin(2),fun(xmin,fnum),'r*');
end

% Contour and quiver plots
dfdx1 = zeros(size(x1));
dfdx2 = zeros(size(x1));
for i = 1 : size(x1,1)
	for j = 1 : size(x1,2)
		x = [x1(i,j),x2(i,j)]';
		[~,J] = fun(x,fnum);
		dfdx1(i,j) = J(1);
		dfdx2(i,j) = J(2);
	end
end
subplot(1,2,2);
contour(x1,x2,f);
hold on;
quiver(x1,x2,dfdx1,dfdx2);
if drawMin
	plot(xmin(1),xmin(2),'r*');
end
axis equal;
axis([min(min(x1)),max(max(x1)),min(min(x2)),max(max(x2))]);
xlabel('x_1');
ylabel('x_2');

% Show where the function is convex, concave, etc.
if drawAreas
	c = zeros(size(x1));
	for i = 1 : size(x1,1)
		for j = 1 : size(x1,2)
			x = [x1(i,j),x2(i,j)]';
			[~,~,H] = fun(x);
			E = eig(H);
			if E(1) >= 0
				if E(2) >= 0
					c(i,j) = 1;
				else
					c(i,j) = 2;
				end
			else
				if E(2) >= 0
					c(i,j) = 3;
				else
					c(i,j) = 4;
				end
			end
		end
	end
	subplot(1,2,1);
	hold off;
	surf(x1,x2,f,c);
	hold on;
	grid on;
	axis equal;
	xlabel('x_1');
	ylabel('x_2');
	zlabel('f');
end
end
