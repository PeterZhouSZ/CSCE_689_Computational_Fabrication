function driver(prob)

if nargin < 1
	prob = 3;
end

iterMax = 100;
thresh = 1e-6;

switch prob
	case 0
		fnum = 0;
		xmin = -1;
		xmax = 3;
		ymin = -4;
		ymax = 0;
		dx = 0.1;
		dy = 0.1;
		momentum = false;
		newton = true;
		xInit = [0 0]';
		alpha = 1;
		beta = 0.2;
		gamma = 1.0;
	case 1
		fnum = 1;
		xmin = -5;
		xmax = 7;
		ymin = -5;
		ymax = 7;
		dx = 0.5;
		dy = 0.5;
		momentum = false;
		newton = true;
		xInit = [5 7]';
		alpha = 1;
		beta = 0.6;
		gamma = 1.0;
	case 2
		fnum = 1;
		xmin = -5;
		xmax = 7;
		ymin = -5;
		ymax = 7;
		dx = 0.5;
		dy = 0.5;
		momentum = false;
		newton = true;
		xInit = [-0.5 -1]';
		alpha = 1.0;
		beta =0.1;
		gamma = 1.0;
    case 3
        fnum = 2;
		xmin = -5;
		xmax = 7;
		ymin = -5;
		ymax = 7;
		dx = 0.5;
		dy = 0.5;
		momentum = false;
		newton = true;
		xInit = [-0.5 -1]';
		alpha = 1;
		beta =0.1;
		gamma = 1.0;
     case 4
        fnum = 2;
		xmin = -5;
		xmax = 7;
		ymin = -5;
		ymax = 7;
		dx = 0.5;
		dy = 0.5;
		momentum = false;
		newton = true;
		xInit = [0 -2.5]';
		alpha = 1;
		beta =0.1;
		gamma = 1.0;
end
draw(fnum,xmin,xmax,dx,ymin,ymax,dy);

% Gradient Descent
  x = xInit;
  f = fun(x,fnum);
  subplot(1,2,1);
  plot3(x(1),x(2),f,'ro:','LineWidth',2,'MarkerSize',5);
  subplot(1,2,2);
  plot(x(1),x(2),'ro:','LineWidth',2,'MarkerSize',5);
for iter = 1 : iterMax
	%Your code here
    [f, g, H] = fun(x,fnum);
    dx = -alpha* g;
    if (norm(dx)<thresh)
        subplot(1,2,1);
        plot3(x(1),x(2),f,'ro:','LineWidth',2,'MarkerSize',5);
        subplot(1,2,2);
        plot(x(1),x(2),'ro:','LineWidth',2,'MarkerSize',5);
        break;
    end
    x = x+dx;  
    subplot(1,2,1);
    plot3(x(1),x(2),f,'ro:','LineWidth',2,'MarkerSize',5);
    subplot(1,2,2);
    plot(x(1),x(2),'ro:','LineWidth',2,'MarkerSize',5);
end
fprintf('Gradient Descent: %d iters\n',iter);

% Gradient Descent with Momentum
if momentum
	v = [0 0]';
	x = xInit;
	for iter = 1 : iterMax
		% Your code here
        [f, g, H] = fun(x,fnum);
        v = beta * v - alpha * g;
        dx = v;
        if (norm(dx) < thresh)
            break;
        end
        x = x + dx;  
        subplot(1,2,1);
        plot3(x(1),x(2),f,'bo:','LineWidth',2,'MarkerSize',5);
        subplot(1,2,2);
        plot(x(1),x(2),'bo:','LineWidth',2,'MarkerSize',5);
	end
	fprintf('Momentum: %d iters\n',iter);
end

% Newton's method
if newton
	x = xInit;
	for iter = 1 : iterMax
		% Your code here
        [f, g, H] = fun(x,fnum);
        dx = - inv(H)*g;
        if (norm(dx) < thresh)
            break;
        end
        x = x + dx;
       [f, g, H] = fun(x,fnum);
        subplot(1,2,1);
        plot3(x(1),x(2),f,'go:','LineWidth',2,'MarkerSize',5);
        subplot(1,2,2);
        plot(x(1),x(2),'go:','LineWidth',2,'MarkerSize',5);
        
	end
	fprintf('Newton: %d iters\n',iter);
end

end
