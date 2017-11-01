function linkages(scene)
if nargin < 1
	scene = 0;
end

links = [];
pins = [];
sliders = [];
particles = [];

% lsqnonlin options
opt = optimoptions('lsqnonlin');
opt.Jacobian = 'on';
opt.DerivativeCheck = 'off';
opt.TolFun = 1e-12;
opt.Display = 'off';

% Set up the scene here.
% Note that links don't have to be placed exactly. The first call to
% solveLinkage() will automatically snap the links so that all the
% constraints are satisfied.

% These values can be overridden inside the switch statement
T = 3; % final time
dt = 0.01; % time step
drawHz = 30; % refresh rate
angVel = 2*pi; % driver angular velocity
oscillate = 0; % specify a range that the driver oscillates between
%oscillate = 0; % integrates the target angle with a constant angular velocity
if oscillate == 0
    angVel = 2*pi;
end
if oscillate == 1
    angle1 = pi/7;
    angle2 = pi/3;
    angVel = 5;
end

switch scene
	case 0
		% Crank-rocker
		% Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Left link
		links(2).angle = pi/2;
		links(2).pos = [-10 0]';
		links(2).verts = [
			-1.0  11.0  11.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Right link
		links(3).angle = pi/2;
		links(3).pos = [10 0]';
		links(3).verts = [
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Top link (note: links don't need to be rectangular)
		links(4).angle = 0;
		links(4).pos = [-10 10]';
		links(4).verts = [
			-1.0  25.0 31.0  31.0  5.0 -1.0
			-1.0  -5.0 -1.0   1.0  3.0  1.0
			];
		
		% Which link is grounded?
		grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
		
		% Bottom-left
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
		% Bottom-right
		pins(end+1).linkA = 1;
		pins(end).linkB = 3;
		pins(end).pointA = [20,0]';
		pins(end).pointB = [0,0]';
		% Left-top
		pins(end+1).linkA = 2;
		pins(end).linkB = 4;
		pins(end).pointA = [10,0]';
		pins(end).pointB = [0,0]';
		% Right-top
		pins(end+1).linkA = 3;
		pins(end).linkB = 4;
		pins(end).pointA = [10+10*rand(1),0]'; % pin location on link3 is randomized
		pins(end).pointB = [20,0]';
		
		% List of tracer particles for display
		particles(1).link = 4; % which link?
		particles(1).point = [5,3]'; % tracer particle point in local coords
		particles(2).link = 4;
		particles(2).point = [25,-5]';
	case 1
        
        % Drag-link
		% Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Left link
		links(2).angle = pi/2;
		links(2).pos = [-10 0]';
		links(2).verts = [
			-1.0  15.0  15.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Right link
		links(3).angle = pi/2;
		links(3).pos = [10 0]';
		links(3).verts = [
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Top link (note: links don't need to be rectangular)
		links(4).angle = 0;
		links(4).pos = [-10 10]';
		links(4).verts = [
			-1.0  25.0 31.0  31.0  5.0 -1.0
			-1.0  -5.0 -1.0   1.0  3.0  1.0
			];
		
		% Which link is grounded?
		grounded = 2;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 4;
		
		% Bottom-left
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
		% Bottom-right
		pins(end+1).linkA = 1;
		pins(end).linkB = 3;
		pins(end).pointA = [20,0]';
		pins(end).pointB = [0,0]';
		% Left-top
		pins(end+1).linkA = 2;
		pins(end).linkB = 4;
		pins(end).pointA = [15,0]';
		pins(end).pointB = [0,0]';
		% Right-top
		pins(end+1).linkA = 3;
		pins(end).linkB = 4;
		pins(end).pointA = [10+10*rand(1),0]'; % pin location on link3 is randomized
		pins(end).pointB = [20,0]';
		
		% List of tracer particles for display
		particles(1).link = 1; % which link?
		particles(1).point = [5,0]'; % tracer particle point in local coords
		particles(2).link = 4;
		particles(2).point = [25,-5]';
	case 2
		% Double-rocker
        oscillate = 1;
        angle1 = pi/7;
        angle2 = pi/3;
        angVel = 2;

        % Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  11.0  11.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Left link
		links(2).angle = pi/2;
		links(2).pos = [-10 0]';
		links(2).verts = [
			-1.0  10.0  10.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Right link
		links(3).angle = pi/2;
		links(3).pos = [0 0]';
		links(3).verts = [
			-1.0  9.0  9.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Top link (note: links don't need to be rectangular)
		links(4).angle = 0;
		links(4).pos = [-10 8]';
		links(4).verts = [
			-1.0  4.0 6.0  6.0  -1.0
			-1.0  -1.0 -1.0   1.0  1.0
			];
		
		% Which link is grounded?
		grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
		
		% Bottom-left
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
		% Bottom-right
		pins(end+1).linkA = 1;
		pins(end).linkB = 3;
		pins(end).pointA = [10,0]';
		pins(end).pointB = [0,0]';
		% Left-top
		pins(end+1).linkA = 2;
		pins(end).linkB = 4;
		pins(end).pointA = [8,0]';
		pins(end).pointB = [0,0]';
		% Right-top
		pins(end+1).linkA = 3;
		pins(end).linkB = 4;
		pins(end).pointA = [8,0]'; % pin location on link3 is randomized
		pins(end).pointB = [4,0]';
		
		% List of tracer particles for display
		particles(1).link = 4; % which link?
		particles(1).point = [0,0]'; % tracer particle point in local coords
		particles(2).link = 4;
		particles(2).point = [4,0]';
		
	case 3
		% Hoekens
		% Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  11.0  11.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Left link
		links(2).angle = pi/2;
		links(2).pos = [-10 0]';
		links(2).verts = [
			-1.0  5.0 5.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Right link
		links(3).angle = pi/2;
		links(3).pos = [10 0]';
		links(3).verts = [
			-1.0  11.0  11.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% Top link (note: links don't need to be rectangular)
		links(4).angle = 0;
		links(4).pos = [-10 10]';
		links(4).verts = [
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0  1.0  1.0
			];
		
		% Which link is grounded?
		grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
		
		% Bottom-left
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
		% Bottom-right
		pins(end+1).linkA = 1;
		pins(end).linkB = 3;
		pins(end).pointA = [8,0]';
		pins(end).pointB = [0,0]';
		% Left-top
		pins(end+1).linkA = 2;
		pins(end).linkB = 4;
		pins(end).pointA = [4,0]';
		pins(end).pointB = [0,0]';
		% Right-top
		pins(end+1).linkA = 3;
		pins(end).linkB = 4;
		pins(end).pointA = [10,0]'; % pin location on link3 is randomized
		pins(end).pointB = [10,0]';
		
		% List of tracer particles for display
		particles(1).link = 4; % which link?
		particles(1).point = [0,0]'; % tracer particle point in local coords
		particles(2).link = 4;
		particles(2).point = [20,0]';
        
        
        
        
	case 4
		% Peaucellier-Lipkin
	case 5
		% Klann
        oscillate = 0;
        
        % Triangle
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [0 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			0.0     26.616  0.0
			6.145   0.0     -13.0
			];
		% crank
		links(2).angle = 3/4*pi;
		links(2).pos = [26.616 0]';
		links(2).verts = [
			-1.0  12.0 12.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% 
		links(3).angle = -pi*4/5;
		links(3).pos = [17 0]';
		links(3).verts = [
			-1.0  29  29  -1.0
			-1.0  -1.0   1.0   1.0
			];
		% 
		links(4).angle = (160/180)*pi;
		links(4).pos = [0 6.145]';
		links(4).verts = [
			-1.0  19.0  19.0 -1.0
			-1.0  -1.0  1.0  1.0
			];
		%
        links(5).angle = -pi*(170/180);
		links(5).pos = [0 -13.0]';
		links(5).verts = [
			-1.0  14.0  14.0 -1.0
			-1.0  -1.0  1.0  1.0
			];
		%
        links(6).angle =-pi*(130/180);
		links(6).pos = [0 -13.0]';
		links(6).verts = [
			-1.0  27.5  27.5 -1.0
			-1.0  -1.0  1.0  1.0
			];
        
        %
        links(7).angle = -pi/5*6;
		links(7).pos = [-13/2*sqrt(2) 13/2*sqrt(2)-13.0]';
		links(7).verts = [
			-1.0  23.5  23.5 -1.0
			-1.0  -1.0  1.0  1.0
			];
        %
        links(8).angle = -pi*(100/180);
		links(8).pos = [-30 -10.0]';
		links(8).verts = [
			-1.0  50  50 -1.0
			-1.0  -1.0  1.0  1.0
			];
        
        links(9).angle = -pi;
		links(9).pos =  [17 0]';
		links(9).verts = [
			-1.0  51.80945  51.80945 -1.0
			-1.0  -1.0  1.0  1.0
			];
        
        links(10).angle = -pi;
		links(10).pos =  [-30 8.0]';
		links(10).verts = [
			-1.0 76  76 -1.0
			-1.0  -1.0  1.0  1.0
			];
        
        
		% Which link is grounded?
		grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
		
		% 1-2
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [26.616,0]';
		pins(end).pointB = [0,0]';
		% 1-4
		pins(end+1).linkA = 1;
		pins(end).linkB = 4;
		pins(end).pointA = [0,6.145]';
		pins(end).pointB = [0,0]';
		% 1-5
		pins(end+1).linkA = 1;
		pins(end).linkB = 5;
		pins(end).pointA = [0,-13.0]';
		pins(end).pointB = [0,0]';
		% 5-7
		pins(end+1).linkA = 5;
		pins(end).linkB = 7;
		pins(end).pointA = [13.0,0]'; 
		pins(end).pointB = [0,0]';
		
        % 2-3
		pins(end+1).linkA = 2;
		pins(end).linkB = 3;
		pins(end).pointA = [11,0]';
		pins(end).pointB = [0,0]';
        
        % 3-7
		pins(end+1).linkA = 3;
		pins(end).linkB = 7;
		pins(end).pointA = [28.8,0]'; 
		pins(end).pointB = [0,0]';
        
        % 7-8
		pins(end+1).linkA = 7;
		pins(end).linkB = 8;
		pins(end).pointA = [22.2,0]';
		pins(end).pointB = [0,0]';
        
         % 4-6
		pins(end+1).linkA = 4;
		pins(end).linkB = 6;
		pins(end).pointA = [18.2,0]';
		pins(end).pointB = [0,0]';
        
        % 6-8
		pins(end+1).linkA = 6;
		pins(end).linkB = 8;
		pins(end).pointA = [26.5,0]'; % pin location on link3 is randomized
		pins(end).pointB = [0,0]';
        
         % 3-9
		pins(end+1).linkA = 3;
		pins(end).linkB = 9;
		pins(end).pointA = [0,0]'; % pin location on link3 is randomized
		pins(end).pointB = [0,0]';
        
         % 7-9
		pins(end+1).linkA = 7;
		pins(end).linkB = 9;
		pins(end).pointA = [22.2,0]'; % pin location on link3 is randomized
		pins(end).pointB = [50.80945,0]';
        
        % 6-10
		pins(end+1).linkA = 6;
		pins(end).linkB = 10;
		pins(end).pointA = [0,0]'; % pin location on link3 is randomized
		pins(end).pointB = [0,0]';
        
         % 8-10
		pins(end+1).linkA = 8;
		pins(end).linkB = 10;
		pins(end).pointA = [49,0]'; % pin location on link3 is randomized
		pins(end).pointB = [75.23829,0]';
        
		% 
		particles(1).link = 8; % which link?
		particles(1).point = [49,0]'; % tracer particle point in local coords
		particles(2).link = 8;
		particles(2).point = [0,0]';
        
        
        
	case 6
		% Simple slider
        opt.Jacobian = 'off';
        oscillate = 0; 
        % Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  21.0  21.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% 
		links(2).angle = pi/2;
		links(2).pos = [10 0]';
		links(2).verts = [
			-1.0  30.0  30.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		%
		links(3).angle = pi/2;
		links(3).pos = [-10 0]';
		links(3).verts = [
			-1.0  11.0  11.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		
		
		% Which link is grounded?
		grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 3;
		
		
		% 2-3
		pins(end+1).linkA = 2;
		pins(end).linkB = 3;
		pins(end).pointA = [30,0]';
		pins(end).pointB = [10,0]';
		
        % 1-3
		pins(end+1).linkA = 1;
		pins(end).linkB = 3;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
        
        sliders(end+1).linkA = 1;
        sliders(end).linkB = 2;
        sliders(end).pointA = [20,0]';
        sliders(end).pointB1 = [0,0]';
        sliders(end).pointB2 = [20,0]';
        
		
		% List of tracer particles for display
		particles(1).link = 2; % which link?
		particles(1).point = [0,-1]'; % tracer particle point in local coords
		particles(2).link = 2;
		particles(2).point = [30,-1]';
        
        
        
	case 10
		% Extra credit!
        % Scissor
        oscillate =1; 
		opt.Jacobian = 'off';
        angle1 = 0;
        angle2 = pi/4;
        angVel = 5;
        % Bottom link
		links(1).angle = 0; % rotation from the positive x-axis
		links(1).pos = [-10 0]'; % position of the center of rotation
		links(1).verts = [ % display vertices
			-1.0  31.0  31.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
		% 
		links(2).angle = pi/4;
		links(2).pos = [-10 0]';
		links(2).verts = [
			-1.0  21*sqrt(2)  21*sqrt(2) -1.0
			-1.0  -1.0   1.0  1.0
			];
		%
		links(3).angle = -pi/4*3;
		links(3).pos = [10 0]';
		links(3).verts = [
			-1.0  21*sqrt(2)  21*sqrt(2) -1.0
			-1.0  -1.0   1.0  1.0
			];
		
        links(4).angle = pi/4;
		links(4).pos = [-10 20]';
		links(4).verts = [
			-1.0  21*sqrt(2)  21*sqrt(2) -1.0
			-1.0  -1.0   1.0  1.0
			];
       
        links(5).angle = -pi/4*3;
		links(5).pos = [10 20]';
		links(5).verts = [
			-1.0  21*sqrt(2)  21*sqrt(2) -1.0
			-1.0  -1.0   1.0  1.0
			];
       
        links(6).angle = 0;
		links(6).pos = [-10 40]';
		links(6).verts = [
			-1.0  31.0  31.0 -1.0
			-1.0  -1.0   1.0  1.0
			];
        
		grounded = 1;
		driver = 2;
		
		% 2-3
		pins(end+1).linkA = 2;
		pins(end).linkB = 3;
		pins(end).pointA = [10*sqrt(2),0]';
		pins(end).pointB = [10*sqrt(2),0]';
		
        % 1-2
		pins(end+1).linkA = 1;
		pins(end).linkB = 2;
		pins(end).pointA = [0,0]';
		pins(end).pointB = [0,0]';
        
         % 3-4
		pins(end+1).linkA = 3;
		pins(end).linkB = 4;
		pins(end).pointA = [20*sqrt(2),0]';
		pins(end).pointB = [0,0]';
        
         % 2-5
		pins(end+1).linkA = 2;
		pins(end).linkB = 5;
		pins(end).pointA = [20*sqrt(2),0]';
		pins(end).pointB = [0,0]';
        
         % 4-5
		pins(end+1).linkA = 4;
		pins(end).linkB = 5;
		pins(end).pointA = [10*sqrt(2),0]';
		pins(end).pointB = [10*sqrt(2),0]';
        
         % 5-6
		pins(end+1).linkA = 5;
		pins(end).linkB = 6;
		pins(end).pointA = [20*sqrt(2),0]';
		pins(end).pointB = [0,0]';
        
        % s1-3
        sliders(end+1).linkA = 3;
        sliders(end).linkB = 1;
        sliders(end).pointA = [0,0]';
        sliders(end).pointB1 = [20,0]';
        sliders(end).pointB2 = [20*sqrt(2),0]';       
		
        % s1-3
        sliders(end+1).linkA = 4;
        sliders(end).linkB = 6;
        sliders(end).pointA = [20*sqrt(2),0]';
        sliders(end).pointB1 = [20,0]';
        sliders(end).pointB2 = [20*sqrt(2),0]';  
        
        particles(1).link = 3; % which link?
		particles(1).point = [20*sqrt(2),0]'; % tracer particle point in local coords
		particles(2).link = 2;
		particles(2).point = [20*sqrt(2),0]';
        
end

% Initialize
for i = 1 : length(links)
	links(i).grounded = (i == grounded); %#ok<*AGROW>
	links(i).driver = (i == driver);
	% These target quantities are only used for grounded and driver links
	links(i).angleTarget = links(i).angle;
	links(i).posTarget = links(i).pos;
end
for i = 1 : length(particles)
	particles(i).pointsWorld = zeros(2,0); % transformed points, initially empty
end

% Debug: drawing here to debug scene setup
%drawScene(0,links,pins,sliders,particles);

% Simulation loop
t = 0; % current time
t0 = -inf; % last draw time

while t < T
	% Procedurally set the driver angle.
	% Right now, the target angle is being linearly increased, but you may
	% want to do something else.
	if oscillate == 0
        links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
    end
    
    if oscillate == 1
        
        ang = t * angVel;
        links(driver).angleTarget = (sin(ang))^2 * angle1 + (1-sin(ang)^2) * angle2;
    end
    
	% Solve for linkage orientations and positions
	[links,feasible] = solveLinkage(links,pins,sliders,opt);
	% Update particle positions
	particles = updateParticles(links,particles);
	% Draw scene
	if t - t0 > 1 / drawHz
		drawScene(t,links,pins,sliders,particles);
		t0 = t;
	end
	% Quit if over-constrained
	if ~feasible
		break;
	end
	t = t + dt;
end
drawScene(t,links,pins,sliders,particles);

end

%%
function [R,dR] = rotationMatrix(angle)
c = cos(angle);
s = sin(angle);
% Rotation matrix
R = zeros(2);
R(1,1) = c;
R(1,2) = -s;
R(2,1) = s;
R(2,2) = c;
if nargout >= 2
	% Rotation matrix derivative
	dR = zeros(2);
	dR(1,1) = -s;
	dR(1,2) = -c;
	dR(2,1) = c;
	dR(2,2) = -s;
end
end

%%
function [links,feasible] = solveLinkage(links,pins,sliders,opt)
nlinks = length(links);
% Extract the current angles and positions into a vector
angPos0 = zeros(3*nlinks,1);
for i = 1 : nlinks
	link = links(i);
	ii = (i-1)*3+(1:3);
	angPos0(ii(1)) = link.angle;
	angPos0(ii(2:3)) = link.pos;
end
% Limits
lb = -inf(size(angPos0));
ub =  inf(size(angPos0));
% Solve for angles and positions
[angPos,r2] = lsqnonlin(@(angPos)objFun(angPos,links,sliders,pins),angPos0,lb,ub,opt);
% If the mechanism is feasible, then the residual should be zero
feasible = true;
if r2 > 1e-6
	fprintf('Mechanism is over constrained!\n');
	feasible = false;
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
	ii = (i-1)*3+(1:3);
	links(i).angle = angPos(ii(1));
	links(i).pos = angPos(ii(2:3));
end

% Check for underconstrained configuration.
% This check only works when there are no sliders.
if isempty(sliders)
	[c,J] = objFun(angPos,links,sliders,pins); %#ok<ASGLU>
	if rank(J,1e-2) < size(J,2)
		fprintf('Mechanism is under-constrained!\n');
	end
end
end

%%
function [c,J] = objFun(angPos,links,sliders,pins)

nlinks = length(links);
npins = length(pins);
nsliders = length(sliders);

% There cannot be any constraints between kinematic links
for i = 1 : length(pins)
	pin = pins(i);
	indLinkA = pin.linkA; % array index of link A
	indLinkB = pin.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	kinematicA = linkA.grounded || linkA.driver;
	kinematicB = linkB.grounded || linkB.driver;
	if kinematicA && kinematicB
		npins = npins - 1;
	end
end
for i = 1 : length(sliders)
	slider = sliders(i);
	indLinkA = slider.linkA; % array index of link A
	indLinkB = slider.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	kinematicA = linkA.grounded || linkA.driver;
	kinematicB = linkB.grounded || linkB.driver;
	if kinematicA && kinematicB
		nsliders = nsliders - 1;
	end
end

% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
	ii = (i-1)*3+(1:3);
	links(i).angle = angPos(ii(1));
	links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins + nsliders; % 3 for ground, 3 for driver, 2*npins, plus 1*nsliders
c = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
	link = links(i);
	if link.grounded || link.driver
		% Grounded and driver links have their angles and positions
		% prescribed.
		c(k+1    ) = link.angle - link.angleTarget;
		c(k+(2:3)) = link.pos - link.posTarget;
		% The Jacobian of this constraint is the identity matrix
		colAng = (i-1)*3 + 1;
		colPos = (i-1)*3 + (2:3);
		J(k+1,    colAng) = 1;
		J(k+(2:3),colPos) = eye(2);
		k = k + 3;
	end
end
% Pin constraints
for i = 1 : length(pins)
	pin = pins(i);
	indLinkA = pin.linkA; % array index of link A
	indLinkB = pin.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	kinematicA = linkA.grounded || linkA.driver;
	kinematicB = linkB.grounded || linkB.driver;
	if kinematicA && kinematicB
		continue;
	end
	rows = k+(1:2); % row index of this pin constraint
	k = k + 2;
	[Ra,dRa] = rotationMatrix(linkA.angle);
	[Rb,dRb] = rotationMatrix(linkB.angle);
	% Local positions
	ra = pin.pointA;
	rb = pin.pointB;
	% World positions
	xa = Ra * ra + linkA.pos;
	xb = Rb * rb + linkB.pos;
	p = xa(1:2) - xb(1:2);
	c(rows) = p;
	%
	% Optional Jacobian computation
	%
	% Column indices for the angles and positions of links A and B
	colAngA = (indLinkA-1)*3 + 1;
	colPosA = (indLinkA-1)*3 + (2:3);
	colAngB = (indLinkB-1)*3 + 1;
	colPosB = (indLinkB-1)*3 + (2:3);
	% The Jacobian of this constraint is the partial derivative of f wrt
	% the angles and positions of the two links.
	J(rows,colAngA) = dRa * ra;
	J(rows,colPosA) = eye(2);
	J(rows,colAngB) = -dRb * rb;
	J(rows,colPosB) = -eye(2);
end
% Slider constraints
for i = 1 : length(sliders)
	% IMPLEMENT ME
    
    slider = sliders(i);
	indLinkA = slider.linkA; % array index of link A
	indLinkB = slider.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	kinematicA = linkA.grounded || linkA.driver;
	kinematicB = linkB.grounded || linkB.driver;
	if kinematicA && kinematicB
		continue;
	end
	rows = k+1; % row index of this slider constraint
	k = k + 1;
	[Ra] = rotationMatrix(linkA.angle);
	[Rb] = rotationMatrix(linkB.angle);
    % Local positions
	ra = slider.pointA;
	rb1 = slider.pointB1;
    rb2 = slider.pointB2;
    % World positions
	xa = Ra * ra + linkA.pos;
	xb1 = Rb * rb1 + linkB.pos;
    xb2 = Rb * rb2 + linkB.pos;
	a = 1/2*((xb1(1)-xa(1))*(xb2(2)-xa(2)) - (xb2(1)-xa(1))*(xb1(2)-xa(2)));
	c(rows) = a;
end
end

%%
function particles = updateParticles(links,particles)
% Transform particle position from local to world
for i = 1 : length(particles)
	particle = particles(i);
	link = links(particle.link);
	% IMPLEMENT ME: compute x, the world space position of the particle.
    x =  rotationMatrix(link.angle)* particle.point + link.pos;
	% Append world position to the array (grows indefinitely)
	particles(i).pointsWorld(:,end+1) = x;
end
end

%%
function drawScene(t,links,pins,sliders,particles)
if t == 0
	clf;
	axis equal;
	hold on;
	grid on;
	xlabel('X');
	ylabel('Y');
end
cla;
% Draw links
rlines = [];
glines = [];
blines = [];
nan2 = nan(2,1);
for i = 1 : length(links)
	link = links(i);
	% Draw frame
	R = rotationMatrix(link.angle);
	p = link.pos; % frame origin
	s = 1.0; % frame display size
	px = p + s*R(:,1); % frame x-axis
	py = p + s*R(:,2); % frame y-axis
	rlines = [rlines, p, px, nan2];
	glines = [glines, p, py, nan2];
	% Draw link geometry
	E = [R,link.pos;0,0,1]; % transformation matrix
	vertsLocal = [link.verts;ones(1,size(link.verts,2))];
	vertsWorld = E * vertsLocal;
	if link.grounded
		rlines = [rlines, vertsWorld(1:2,[1:end,1]), nan2];
	elseif link.driver
		glines = [glines, vertsWorld(1:2,[1:end,1]), nan2];
	else
		blines = [blines, vertsWorld(1:2,[1:end,1]), nan2];
	end
end
plot(rlines(1,:),rlines(2,:),'r');
plot(glines(1,:),glines(2,:),'g');
plot(blines(1,:),blines(2,:),'b');
% Draw pins
xas = [];
xbs = [];
for i = 1 : length(pins)
	pin = pins(i);
	linkA = links(pin.linkA);
	linkB = links(pin.linkB);
	Ra = rotationMatrix(linkA.angle);
	Rb = rotationMatrix(linkB.angle);
	xa = Ra * pin.pointA + linkA.pos;
	xb = Rb * pin.pointB + linkB.pos;
	xas = [xas, xa];
	xbs = [xbs, xb];
end
plot(xas(1,:),xas(2,:),'co','MarkerSize',10,'MarkerFaceColor','c');
plot(xbs(1,:),xbs(2,:),'mx','MarkerSize',10,'LineWidth',2);
% Draw sliders
xaas = [];
xxs = [];
yys =[];
if ~isempty(sliders)
	% IMPLEMENT ME
    for i = 1 : length(sliders)
        slider = sliders(i);
        linkA = links(slider.linkA);
        linkB = links(slider.linkB);
        Ra = rotationMatrix(linkA.angle);
        Rb = rotationMatrix(linkB.angle);
        xa = Ra * slider.pointA + linkA.pos;
        xb1 = Rb * slider.pointB1 + linkB.pos;
        xb2 = Rb * slider.pointB2 + linkB.pos;
        xaas = [xaas, xa];
        xxs = [xxs, xb1(1)];
        xxs = [xxs, xb2(1)];
        yys = [yys, xb1(2)];
        yys = [yys, xb2(2)];
        plot(xxs,yys,'m','LineWidth',3);
        xxs = [];
        yys =[];
    end
    plot(xaas(1,:),xaas(2,:),'co','MarkerSize',10,'MarkerFaceColor','c');
end


% Draw particles
ps = [];
ps1 = [];
for i = 1 : length(particles)
	particle = particles(i);
	if ~isempty(particle.pointsWorld)
		ps = [ps, particle.pointsWorld];
		ps1 = [ps1, particle.pointsWorld(:,end)];
	end
	ps = [ps, nan2];
end
plot(ps(1,:),ps(2,:),'k');
plot(ps1(1,:),ps1(2,:),'ko');
% End draw
title(sprintf('t=%.3f',t));
drawnow;
end
