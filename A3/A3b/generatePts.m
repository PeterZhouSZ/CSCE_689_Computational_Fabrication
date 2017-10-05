function points = generatePts(caseNumber,count,args)

% Create points
points = zeros(2,count);
for i = 1 : count
	% s is in [0,1]
	s = (i-1.0)/(count-1.0);
	switch caseNumber
		case 0
			p = parametricFunction0(s,args);
		case 1
			p = parametricFunction1(s,args);
        case 2
            p = parametricFunction2(s,args);
            
	end
	points(:,i) = p;
end

% Subtract center of mass
points = points - mean(points,2);

end

function p = parametricFunction0(s,args)
A = args(1);
a = args(2);
x = s;
y = A*sin(a*(2.0*pi)*s);
p = [x,y]';
end

function p = parametricFunction1(s,args)
% https://en.wikipedia.org/wiki/Lissajous_curve
A = args(1);
B = args(2);
a = args(3);
b = args(4);
d = args(5);
t = 2.0*pi*s;
x = A*sin(a*t+d);
y = B*sin(b*t);
p = [x,y]';
end

function p = parametricFunction2(s,args)
% Maurer Rose
%
a = args(1);
n = args(2);
d = args(3);
theta = 360*s;
r = a*sin(n*theta*d);
x = r*cos(theta);
y = r*sin(theta);
p = [x,y]';
end