function driver(caseNumber)

if nargin < 1
	caseNumber = 2;
end

% Use 'shuffle' if you want a different randomized input every time.
% Use 0 if you want the same randomized input every time.
rng(0);
%rng(0);

% Default perturbation parameters
% These can be changed inside the switch statement if necessary.
%theta = randBetween(-50.0,50.0)*pi/180.0;
theta =(70)*pi/180;
%theta = randBetween(70.0,120.0)*pi/180.0;

x = randBetween(-1.0,1.0);
y = randBetween(-1.0,1.0);
noise = 0.01;

% Create source and target data points
switch caseNumber
	case 0
		targetPts = generatePts(caseNumber,100,[0.2 2.0]);
		sourcePts = generatePts(caseNumber,200,[0.2 2.0]);
	case 1
		targetPts = generatePts(caseNumber,100,[1.0 1.0 3.0 2.0 0.0]);
		sourcePts = generatePts(caseNumber,200,[1.0 1.0 3.0 2.0 0.0]);
    case 2
        %targetPts = generatePts(caseNumber,360,[1 1 2 39]);
        %sourcePts = generatePts(caseNumber,360,[1 1 2 39]);
       
        % 3 petals
        targetPts = generatePts(caseNumber,360,[1 1 3 47]);
        sourcePts = generatePts(caseNumber,360,[1 1 3 47]);
        
        % 4 petals
        %targetPts = generatePts(caseNumber,360,[1 1 2 39]);
        %sourcePts = generatePts(caseNumber,360,[1 1 2 39]);
        
        % 7 petals
        %targetPts = generatePts(caseNumber,360,[1 1 7 19]);
        %sourcePts = generatePts(caseNumber,360,[1 1 7 19]);
        
        % 8 petals
        %targetPts = generatePts(caseNumber,360,[1 1 4 31]);
        %sourcePts = generatePts(caseNumber,360,[1 1 4 31]);
        
        % 12 petals
        %targetPts = generatePts(caseNumber,360,[1 1 6 71]);
        %sourcePts = generatePts(caseNumber,360,[1 1 6 71]);
        
       
end

% Perturb data
sourcePts = perturbPts(sourcePts,theta,x,y,noise);

% Plot source and target
clf
hold on;
plot(sourcePts(1,:),sourcePts(2,:),'ro:','LineWidth',1);
plot(targetPts(1,:),targetPts(2,:),'gx-','LineWidth',2);
axis equal

% Run ICP
thresh = 1e-5;
iterMax = 100;
nsamples = 360;
medianMult = inf;
[thetaXY,iter,total_reject] = icp2d(sourcePts,targetPts,thresh,iterMax,nsamples,medianMult);
total_reject
% Plot results
sourcePts1 = transformPts(thetaXY,sourcePts);
plot(sourcePts1(1,:),sourcePts1(2,:),'bo');
title(sprintf('%d iterations\n',iter));

end

function r = randBetween(r0,r1)
s = rand();
r = (1.0-s)*r0 + s*r1;
end
