function points1 = transformPts(thetaXY,points)
% Transform by current theta, x, and y
sourceCount = size(points,2);
s = sin(thetaXY(1));
c = cos(thetaXY(1));
R = [c,-s;s,c];
t = thetaXY(2:3);
% Vectorized version for efficiency, since for-loops are slow in Matlab
points1 = R*points + repmat(t,1,sourceCount);
end
