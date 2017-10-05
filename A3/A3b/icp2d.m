function [thetaXY,iter,total_reject] = icp2d(sourcePts,targetPts,thresh,iterMax,nsamples,medianMult)
% sourcePts: input source points
% targetPts: input target points
% thresh: matching threshold
% iterMax: maximum number of iterations
% nsamples: number of samples to choose from the source points
% medianMult: the median multiplier for rejecting pairs

% Optimizer options: shouldn't need to modify
opt = optimoptions('fminunc');
opt.SpecifyObjectiveGradient = true;
opt.CheckGradients = false;
opt.FiniteDifferenceType = 'central';
opt.Display = 'off';

% Number of points
sourceCount = size(sourcePts,2);
targetCount = size(targetPts,2);

% Generate sample points
p = randperm(sourceCount, nsamples);
sampleSourcePts = sourcePts(:,p);

% Current theta, x, and y, stored as a column vector:
% thetaXY(1) is theta
% thetaXY(2) is x
% thetaXY(3) is y
thetaXY = [0 0 0]';
closest = zeros(nsamples,1);
thetaXY1 = thetaXY;
distList = zeros(nsamples,1);
total_reject=0;

for iter = 1 : iterMax
	% Do something here
    % Compute the new transformed points using the current values of
    % thetaXY
    thetaXY = thetaXY1;
    sampleSourcePts1 = transformPts(thetaXY, sampleSourcePts);
    nsamples = size(sampleSourcePts,2);
    % Compute the closest targetPts to all sample points
    for i = 1 : nsamples
        dist_min = inf;
        for j = 1 : targetCount
            dist = norm(sampleSourcePts1(:,i)-targetPts(:,j));
            if(dist < dist_min)
                dist_min = dist;
                closest(i) = j;
                distList(i) = dist;
            end
        end
    end
    
    medianDist = median(distList);
    reject = [];
    for i = 1 : nsamples
        if(distList(i) > medianMult * medianDist)
            reject = [reject, i];
        end
    end
    total_reject=total_reject + size(reject,2);
    sampleSourcePts(:,reject)=[];
    closest(reject)=[];
    distList(reject)=[];
    
    thetaXY1 = fminunc(@(thetaXY)distFun(thetaXY,sampleSourcePts,targetPts,closest),thetaXY,opt);

    % Test for convergence
    if(norm(thetaXY1 - thetaXY) < thresh)
        break;
    end
end

end

function [f,J] = distFun(thetaXY,sourcePts,targetPts,closest)
% Transform sourcePts by the current theta, x, and y.
sourcePts = transformPts(thetaXY,sourcePts);
% Find the total distance between the corresponding points in sourcePts1
% and targetPts.
Dx = sourcePts - targetPts(:,closest);
f = 0.5*sum(sum(Dx.^2));
% Gradients:
% df/dTheta = df/ddx*ddx/dxi*dxi/dTheta
%           =   dx' *  (1)  *dR*sourcePt
% dR = [-s,-c;c,-s]
% df/dX = df/ddx*ddx/dxi*dxi/dX
%       =   dx' *  (1)  *[1 0]'
if nargout == 2
	J = zeros(3,1);
	s = sin(thetaXY(1));
	c = cos(thetaXY(1));
	dR = [-s,-c;c,-s];
	J(1) = sum(sum(Dx.*(dR*sourcePts)));
	J(2:3) = sum(Dx,2);
end
end

