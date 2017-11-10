function fem2D(scene)
% FEM explicit triangles
% model = 1; % StVK 

if nargin < 1
	scene = 5;
end

video = true;
if video == true
	video = VideoWriter('output','MPEG-4');
	video.open();
else
	video = [];
end

switch(scene)
	case 0
		nTiles = 8; % Number of tiles in x and y
		dt = 1e-2; % time step
		tEnd = 1.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 5e0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % Neo Hookean
        
	case 1
        nTiles = 1; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % Neo Hookean 
	case 2
        nTiles = 8; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 3.5e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 1; % StVK 
	case 3
        nTiles = 8; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 2.5e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 1; % StVK
	case 4
        nTiles = 8; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 0.3e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % NeoH
	case 5
        nTiles = 8; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 2.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % NeoH
	case 6
        % Unit circle
        n = 1; 
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % NeoH
     case 7
         % scene 5 (b)
        nTiles = 8; % Number of tiles in x and y
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 4.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        E2 = 2e3; % Young's modulus
		nu2 = 0.4; % Poisson's ratio 
        model = 0; % NeoH   
     case 8
        % Polygon
        n = 2; 
		dt = 1e-3; % time step
		tEnd = 3.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 2e0; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
        model = 0; % NeoH
end

% Convert to lambda and mu
% ### TODO ###

lambda = E * nu /((1.0 + nu)*(1.0 - 2.0*nu));
mu = E /(2.0*(1.0 + nu));

% Creates triangles from a regular grid of nodes
if(scene ~= 6 && scene ~= 8)
    [nodes,tris] = createSquare(nTiles);
end
if(scene == 6 || scene == 8)
    [nodes,tris] = createTriangles(n);
end
nNodes = length(nodes);
nTris = length(tris);

% Find some nodes to fix
 if(scene == 1)
     for k = 1 : nNodes
        if nodes(k).X(2) < -0.9
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
 end
 if(scene == 2 || scene == 3 ||scene == 7)
     for k = 1 : nNodes
        if nodes(k).X(2) > 0.9
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
end
 if(scene == 4)
     for k = 1 : nNodes
        if nodes(k).X(1) < -0.9
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
 end
 
 % Fix the middle nodes
if(scene == 5)
     for k = 1 : nNodes
        if nodes(k).X(2) < 0.3 && nodes(k).X(2)>-0.3 && nodes(k).X(1)>-0.3 &&nodes(k).X(1)<0.3
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
end

% Set up material parameters of each triangle
if(scene == 7)
     for k = 1 : nTris
        if k < nTris/4
            tris(k).materials = [E2, nu2];
        else
            tris(k).materials = [E, nu];
        end
    end
end
 
 if(scene == 6)
     for k = 1 : nNodes
        if nodes(k).X(2) > 0.7
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
 end

 if(scene == 8)
     for k = 1 : nNodes
        if nodes(k).X(2) > 0.5
            nodes(k).fixed = true;
        else
            nodes(k).fixed = false;
        end
    end
end
 
% Compute triangle mass and distribute to vertices
% ### TODO ###
for k = 1 : nTris
	tri = tris(k).nodes;
	Xa = nodes(tri(1)).X;
	Xb = nodes(tri(2)).X;
	Xc = nodes(tri(3)).X;
    % Precompute inverse
    Xs{k} = inv([Xb-Xa, Xc-Xa]);
    area = 1/2*((Xb(1)-Xa(1))*(Xc(2)-Xa(2)) - (Xc(1)-Xa(1))*(Xb(2)-Xa(2)));
    triMass = area * rho / 3.0;
	nodes(tri(1)).m = nodes(tri(1)).m + triMass;
	nodes(tri(2)).m = nodes(tri(2)).m + triMass;
	nodes(tri(3)).m = nodes(tri(3)).m + triMass;
end

% Simulation loop
t0 = -inf;
for t = 0 : dt : tEnd
	% Draw scene
	if t - t0 > 1 / drawHz
		draw(t,nodes,tris,video);
		t0 = t;
	end

	% Gravity force
	for k = 1 : nNodes
		nodes(k).f = nodes(k).m*grav;
    end
	
    I = eye(2);
	% FEM force
	% ### TODO ###
    for k = 1 : nTris
        tri = tris(k).nodes;
        xa = nodes(tri(1)).x;
        xb = nodes(tri(2)).x;
        xc = nodes(tri(3)).x;
        dx = [xb - xa, xc - xa];
        
        % Deformation gradient
        F = dx*Xs{k}; 
        
        if(scene == 7)
            E = tris(k).materials(1);
            nu = tris(k).materials(2);
            lambda = E * nu /((1.0 + nu)*(1.0 - 2.0*nu));
            mu = E /(2.0*(1.0 + nu));
        end
        
        % First Piola-Kirchoff stress P
        if(model == 1) 
           % Green strain epsilon
           strain = 1/2 *(F.' * F - I);
           P = F *(2 * mu * strain + lambda * trace(strain) * I);
        end
        
        if(model == 0)
            J = det(F);
            P = mu * (F - inv(F.')) + (lambda * log(J))*inv(F.');   
        end
        % Convert to cauchy stress 
        cauchy = det(inv(F)) * P * inv(F.');
        tris(k).stress = det(inv(F)) * P * inv(F.');
        % Edge normal
        Xa = nodes(tri(1)).X;
        Xb = nodes(tri(2)).X;
        Xc = nodes(tri(3)).X;
        
        ab = Xb - Xa;
        bc = Xc - Xb;
        ca = Xa - Xc;
        
        ab_n = [-ab(2); ab(1)];
        bc_n = [-bc(2); bc(1)];
        ca_n = [-ca(2); ca(1)];
       
        % Normalize 
        ab_n = ab_n / norm(ab_n);
        bc_n = bc_n / norm(bc_n);
        ca_n = ca_n / norm(ca_n);
        
        % Distribute forces to each node
        f_ab = norm(ab) * cauchy * ab_n;
        nodes(tri(1)).f = nodes(tri(1)).f + f_ab * 1/2;
        nodes(tri(2)).f = nodes(tri(2)).f + f_ab * 1/2;
   
        f_bc = norm(bc) * cauchy * bc_n;
        nodes(tri(2)).f = nodes(tri(2)).f + f_bc * 1/2;
        nodes(tri(3)).f = nodes(tri(3)).f + f_bc * 1/2;
    
        f_ca = norm(ca) * cauchy * ca_n;
        nodes(tri(1)).f = nodes(tri(1)).f + f_ca * 1/2;
        nodes(tri(3)).f = nodes(tri(3)).f + f_ca * 1/2;
    end   
	
	% Integrate velocity and position
	% ### TODO ###
    for k = 1 : nNodes
        if(nodes(k).fixed == false)
            nodes(k).v = (nodes(k).m * nodes(k).v + dt * (nodes(k).f))/(nodes(k).m + dt * damping * nodes(k).m);
            nodes(k).x = nodes(k).x + dt * nodes(k).v;
        end
    end
	
end
draw(t,nodes,tris,video);

if ~isempty(video)
	video.close();
end

end

%%
function draw(t,nodes,tris,video)

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
	axis equal;
	%axis(1.5*[-1 1 -1 1]); % Change axis limits here
	grid on;
	view(2);
	colormap jet;
	%caxis([0 50]); % Change color limits here (comment out for auto)
	cb = colorbar;
	ylabel(cb, 'von Mises stress');
end
cla;
hold on;

% Plot triangles with von Mises stress colors
% https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
% (General Plane Stress)
x = [nodes.x]'; % flattened positions
f = reshape([tris.nodes],3,length(tris))'; % flattened indices
stress = reshape([tris.stress],4,length(tris)); % flattened stress
s11 = stress(1,:);
s12 = stress(2,:);
s22 = stress(4,:);
vms = sqrt(s11.*s11 - s11.*s22 + s22.*s22 + 3.0*s12.*s12);
patch('Faces',f,'Vertices',x,'FaceVertexCData',vms','FaceColor','flat');

% Plot fixed points
fixed = [nodes.fixed];
plot(x(fixed,1),x(fixed,2),'ro');

str = sprintf('t = %.4f', t);
title(str);
drawnow;

if ~isempty(video)
	frame = getframe(gcf);
	video.writeVideo(frame);
end

end

%%
function [nodes,tris] = createSquare(n)

% Regular grid with center points
x = linspace(-1,1,n+1);
y = linspace(-1,1,n+1);
[X,Y] = meshgrid(x,y);
x = reshape(X,(n+1)*(n+1),1);
y = reshape(Y,(n+1)*(n+1),1);
dx = 1/n;
dy = 1/n;
xc = linspace(-1+dx,1-dx,n);
yc = linspace(-1+dy,1-dy,n);
[Xc,Yc] = meshgrid(xc,yc);
xc = reshape(Xc,n*n,1);
yc = reshape(Yc,n*n,1);
x = [x;xc];
y = [y;yc];

nodes = [];
for i = 1 : length(x)
	nodes(i).X = [x(i),y(i)]'; %#ok<*AGROW>
	nodes(i).x = nodes(i).X;
	nodes(i).v = [0 0]';
	nodes(i).m = 0;
	nodes(i).f = [0 0]';
	nodes(i).fixed = false;
end

ng = (n+1)*(n+1);
tris = [];
for i = 1 : n
	% Index of the lower left node of the ith column
	ki = (i-1)*(n+1) + 1;
	% Index of the first center node of the ith column
	kci = ng + (i-1)*n + 1;
	for j = 1 : n
		% Index of the lower left node of the jth row of the ith column
		kij = ki + j - 1;
		% Index of the center node of the jth row of the ith column
		kcij = kci + j - 1;
		% Create the four triangles
		tris(end+1).nodes = [kij,kij+(n+1),kcij];
		tris(end+1).nodes = [kij+(n+1),kij+(n+1)+1,kcij];
		tris(end+1).nodes = [kij+(n+1)+1,kij+1,kcij];
		tris(end+1).nodes = [kij+1,kij,kcij];
	end
end
for k = 1 : length(tris)
	tris(k).stress = zeros(2);
end

end

%%
function [nodes,tris] = createTriangles(n)
    % Uniform Mesh on Unit Circle
    if n==1    
    fd=@(p) sqrt(sum(p.^2,2)) -1;
        [p,t] = distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
    end
    % Rectangle with circular hole, refined at circle boundary
    if n==2
       pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
        1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
    [p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
    end
    nodes = [];
    for i = 1 : length(p)
        nodes(i).X = p(i,:)'; 
        nodes(i).x = nodes(i).X;
        nodes(i).v = [0 0]';
        nodes(i).m = 0;
        nodes(i).f = [0 0]';
        nodes(i).fixed = false;
    end

    tris = [];
    for i = 1 : length(t)
        % Create triangles
        tris(end+1).nodes = t(i,:);		
    end
    for k = 1 : length(tris)
        tris(k).stress = zeros(2);
        tris(k).materials = [0 0];
    end

end
