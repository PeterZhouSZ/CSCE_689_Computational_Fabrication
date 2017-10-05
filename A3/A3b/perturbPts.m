function points = perturbPts(points,theta,x,y,noise)
% Rotation and translation
c = cos(theta);
s = sin(theta);
R = [c,-s;s,c];
t = [x;y];
for i = 1 : size(points,2)
	% Apply noise
	points(:,i) = points(:,i) + noise*randn(2,1);
	% Apply rotation and translation
	points(:,i) = R*points(:,i) + t;
end
end
