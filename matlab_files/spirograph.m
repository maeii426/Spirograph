% CSCE 489 Computer Fabrication
% Professor: Dr. Shinjiro Sueda
% Final Project- Fall 2016
% Texas A&M University
%
% Margaret Baxter
% Christine Russell

function spirograph(scene)
workspace;
if nargin < 1
	scene = 1;
end

global video;
video = [];

% Set up the scene with default values
links = [];
pins = [];
sliders = [];
particles = [];
% Default variable values
oscillates = false; 
%multipleRotations = false;
%numR = 0;
clockwise = false;
rangeLower = -pi/2;
rangeUpper = pi/2;
many = false;
num = 1;
secondbase = -1;
switch scene
    case 1
        % Joe Freedman's drawing machine simple config
        % Grounded link
        links(1).grounded = true;
		links(1).angle = 0;
        links(1).pos = [-5 -1]';
        links(1).verts = [
            0.0  9.25 9.25 0.0
           -0.1 -0.1  4.0  4.0
           ];
        %  Driver link
        links(2).grounded = false;
        links(2).angle = -pi;
        links(2).pos = [4.26 0]';
        links(2).verts = [
            0.0  1.0  1.0  0.0
           -0.1 -0.1  0.1  0.1
           ];
        % Arm link
        links(3).grounded = false;
		links(3).angle = -pi; % rotation from the positive x-axis
		links(3).pos = [3.46 0]'; % position of the center of rotation
		links(3).verts = [ % display vertices
			 0.0  10.4 10.4 0.0
			-0.1 -0.1  0.1  0.1
			];
        % Base link
        links(4).grounded = false;
        links(4).angle = 0;
        links(4).pos = [0 0]';
        links(4).verts = [
           0.0  3.6  3.6  0.0
          -0.1 -0.1  0.1  0.1 
          ];
        
        % Which link is grounded?
		% grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
        %secondDriver = 9;
        % Which link is the rotating base?
        base = 4;
        
        % Ground-to-Driver Pin
        pins(1).linkA = 1;
		pins(1).linkB = 2;
		pins(1).pointA = [9.25,1]';
		pins(1).pointB = [0,0]';
        % Driver-to-Arm
        pins(2).linkA = 2;
        pins(2).linkB = 3;
        pins(2).pointA = [0.8 0]';
        pins(2).pointB = [0 0]';
        % Ground-to-Base Pin
        pins(3).linkA = 1;
        pins(3).linkB = 4;
        pins(3).pointA = [5 1]';
        pins(3).pointB = [0 0]';
        
        % Arm Slider pin
        sliders(1).linkA = 3;
        sliders(1).linkB = 1;
        sliders(1).pointA = [10.0 0]';
        sliders(1).pointB = [0.0 0.0]';
        sliders(1).linkAstart = [0.4 0]';
        sliders(1).linkAend = [10.4 0]';
		
        % List of tracer particles for display
		particles(1).link = 3; % which link?
		particles(1).point = [2.8 0.4]'; % tracer particle point in local coords

        many = true;
        num = 14;
    case 2
        % Joe Freedman's drawing machine simple config
        % Grounded link
        links(1).grounded = true;
		links(1).angle = 0;
        links(1).pos = [-5 -1]';
        links(1).verts = [
            0.0  9.25 9.25 0.0
           -0.1 -0.1  4.0  4.0
           ];
        %  Driver link
        links(2).grounded = false;
        links(2).angle = -pi;
        links(2).pos = [4.76 0]';
        links(2).verts = [
            0.0  1.5  1.5  0.0
           -0.1 -0.1  0.1  0.1
           ];
        % Arm link
        links(3).grounded = false;
		links(3).angle = -pi; % rotation from the positive x-axis
		links(3).pos = [3.46 0]'; % position of the center of rotation
		links(3).verts = [ % display vertices
			 0.0  10.4 10.4 0.0
			-0.1 -0.1  0.1  0.1
			];
        % Base link
        links(4).grounded = false;
        links(4).angle = 0;
        links(4).pos = [0 0]';
        links(4).verts = [
           0.0  3.26 3.26 0.0
          -0.1 -0.1  0.1  0.1 
          ];
        
        % Which link is grounded?
		% grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
        %secondDriver = 9;
        % Which link is the rotating base?
        base = 4;
        
        % Ground-to-Driver Pin
        pins(1).linkA = 1;
		pins(1).linkB = 2;
		pins(1).pointA = [9.76,1]';
		pins(1).pointB = [0,0]';
        % Driver-to-Arm
        pins(2).linkA = 2;
        pins(2).linkB = 3;
        pins(2).pointA = [1.3 0]';
        pins(2).pointB = [0 0]';
        % Ground-to-Base Pin
        pins(3).linkA = 1;
        pins(3).linkB = 4;
        pins(3).pointA = [5 1]';
        pins(3).pointB = [0 0]';
        
        % Arm Slider pin
        sliders(1).linkA = 3;
        sliders(1).linkB = 1;
        sliders(1).pointA = [10.0 0]';
        sliders(1).pointB = [1.0 1.0]';
        sliders(1).linkAstart = [0.4 0]';
        sliders(1).linkAend = [10.4 0]';
		
        % List of tracer particles for display
		particles(1).link = 3; % which link?
		particles(1).point = [3.1 -0.4]'; % tracer particle point in local coords

        many = true;
        num = 14;
    case 3
        % Joe Freedman's drawing machine simple config
        % Grounded link
        links(1).grounded = true;
		links(1).angle = 0;
        links(1).pos = [-5 -1]';
        links(1).verts = [
            0.0  9.25 9.25 0.0
           -0.1 -0.1  4.0  4.0
           ];
        %  Driver link
        links(2).grounded = false;
        links(2).angle = -pi;
        links(2).pos = [5 0]';
        links(2).verts = [
            0.0  1.0  1.0  0.0
           -0.1 -0.1  0.1  0.1
           ];
        % Arm link
        links(3).grounded = false;
		links(3).angle = -pi; % rotation from the positive x-axis
		links(3).pos = [4 0]'; % position of the center of rotation
		links(3).verts = [ % display vertices
			 0.0  10.4 10.4 0.0
			-0.1 -0.1  0.1  0.1
			];
        % Base link
        links(4).grounded = false;
        links(4).angle = 0;
        links(4).pos = [0 0]';
        links(4).verts = [
           0.0  4.0  4.0  0.0
          -0.1 -0.1  0.1  0.1 
          ];
        
        % Which link is grounded?
		% grounded = 1;
		% Which link is the driver?
		% Note: the driver must be attached (with a pin) to the ground.
		driver = 2;
        %secondDriver = 9;
        % Which link is the rotating base?
        base = 4;
        
        % Ground-to-Driver Pin
        pins(1).linkA = 1;
		pins(1).linkB = 2;
		pins(1).pointA = [10,1]';
		pins(1).pointB = [0,0]';
        % Driver-to-Arm
        pins(2).linkA = 2;
        pins(2).linkB = 3;
        pins(2).pointA = [0.8 0]';
        pins(2).pointB = [0 0]';
        % Ground-to-Base Pin
        pins(3).linkA = 1;
        pins(3).linkB = 4;
        pins(3).pointA = [5 1]';
        pins(3).pointB = [0 0]';
        
        % Arm Slider pin
        sliders(1).linkA = 3;
        sliders(1).linkB = 1;
        sliders(1).pointA = [10.0 0]';
        sliders(1).pointB = [0.0 0.0]';
        sliders(1).linkAstart = [0.4 0]';
        sliders(1).linkAend = [12.4 0]';
		
        % List of tracer particles for display
		particles(1).link = 3; % which link?
		particles(1).point = [2.8 0.4]'; % tracer particle point in local coords

        many = true;
        num = 14;
end



% Copied some elements from previous homework 4, linkages.m
% Initialize
for i = 1 : length(links)
	%links(i).grounded = (i == grounded); %#ok<*AGROW>
	links(i).driver = (i == driver);
    links(i).base = (i == base);
    links(i).secondbase = (i == secondbase)
    %links(i).secondDriver = (i == secondDriver);
	% These target quantities are only used for grounded and driver links
	links(i).angleTarget = links(i).angle;
	links(i).posTarget = links(i).pos;
end 
for i = 1 : length(particles)
	particles(i).pointsWorld = zeros(2,0); % transformed points, initially empty
    particles(i).pointsBase = zeros(2,0); %relative to the base link
    particles(i).pointsTry = zeros(2,0);
    link = links(particles(i).link);
    R = rotationMatrix(link.angle);
    Rbase = rotationMatrix(links(base).angle);
    x = R * particles(i).point + link.pos;
    xb = inv(Rbase) * x ;
    particles(i).initPos = xb;
end

links(base).pointsA = zeros(2,0);

% Debug: drawing here to debug scene setup
drawScene(0,links,pins,sliders,particles,2*pi,false);

% lsqnonlin options
if verLessThan('matlab','8.1')
    opt = optimset(...
        'Jacobian','off',...
        'DerivativeCheck','off',...
        'Display','off');
else
    opt = optimoptions('lsqnonlin',...
        'Jacobian','off',...
        'DerivativeCheck','off',...
        'Display','off');
end

% Simulation loop
t = 0; % current time
T = 1; % final time
if many
   T = num; 
end
dt = 0.01; % time step
angVel = 2*pi;
%if oscillates %rangeUpper and rangeLower must be specified
%    angVel = (rangeUpper - rangeLower)/2;
%    startingPoint = links(driver).angleTarget + (rangeUpper - angVel);
%else
if clockwise
    angVel = -1*angvel;
end
    
 %   if multipleRotations
 %       angVel = numR*2*pi;
 %   end
%end

while ~(t > T)
	% Procedurally set the driver angle.
	% Right now, the target angle is being linearly increased, but you may
	% want to do something else.
%    if oscillates
%        links(driver).angleTarget = startingPoint + angVel*sin(2*pi*(t+dt));
%    else
     rdriver = links(driver).verts(1,2)-links(driver).verts(2,2)-0.1;
     rbase = links(base).verts(1,2)-links(base).verts(2,2)-0.1;
  %   rsecondbase = links(secondbase).verts(1,2)-links(secondbase).verts(2,2)-0.1;
     links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
     links(base).angleTarget = links(base).angleTarget - (dt*angVel*rdriver)/rbase;
 %    links(secondbase).angleTarget = links(secondbase).angleTarget + (dt*angVel*rdriver)/rsecondbase;
%    end
    
 %   if multipleRotations
 %      links(secondDriver).angleTarget = floor((links(driver).angleTarget)/(2*pi))*pi/(numR*6); 
 %   end
	% Solve for linkage orientations and positions
	[links,feasible] = solveLinkage(links,pins,sliders,opt);
	% Update particle positions
	particles = updateParticles(links,particles,base,t,dt,angVel);
	% Draw scene
	drawScene(t,links,pins,sliders,particles,angVel,false);
	% Quit if over-constrained
	if ~feasible
		break;
	end
	t = t + dt;
end

drawScene(1,links,pins,sliders,particles,2*pi,true)


if ~isempty(video)
	video.close();
end

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
if r2 > 1e-4
	fprintf('Using lenient threshold!\n');
    if r2 > 1e-2
        fprintf('Mechanism is over constrained\n');
        %feasible = false;
    end
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
	ii = (i-1)*3+(1:3);
	links(i).angle = angPos(ii(1));
	links(i).pos = angPos(ii(2:3));
end
end

%%
function [c,J] = objFun(angPos,links,sliders,pins)
nlinks = length(links);
npins = length(pins);
nsliders = length(sliders);
% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
	ii = (i-1)*3+(1:3);
	links(i).angle = angPos(ii(1));
	links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins + length(sliders); % 3 for ground, 3 for driver, 2*npins for pins
c = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
	link = links(i);
	if link.grounded || link.driver || link.base || link.secondbase
		% Grounded and driver links have their angles and positions prescribed.
		c(k+1,    1) = link.angle - link.angleTarget;
		c(k+(2:3),1) = link.pos - link.posTarget;
		k = k + 3;
	end
end
% Pin constraints
for i = 1 : npins
	pin = pins(i);
	rows = k+(1:2); % row index of this pin constraint
	k = k + 2;
	indLinkA = pin.linkA; % array index of link A
	indLinkB = pin.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	[Ra,dRa] = rotationMatrix(linkA.angle);
	[Rb,dRb] = rotationMatrix(linkB.angle);
	% Local positions
	ra = pin.pointA;
	rb = pin.pointB;
	% World positions
	xa = Ra * ra + linkA.pos;
	xb = Rb * rb + linkB.pos;
	p = xa(1:2) - xb(1:2);
	c(rows,1) = p;
end
% Slider constraints
for i = 1 : nsliders
	sl = sliders(i);
	rows = k+(1); % row index of this pin constraint
	k = k + 1;
	indLinkA = sl.linkA; % array index of link A
	indLinkB = sl.linkB; % array index of link B
	linkA = links(indLinkA);
	linkB = links(indLinkB);
	[Ra,dRa] = rotationMatrix(linkA.angle);
	[Rb,dRb] = rotationMatrix(linkB.angle);
	% Local positions
    %ra = sl.pointA;
	rb = sl.pointB;
	% World positions
	%xa = Ra * ra + linkA.pos;
	xb = Rb * rb + linkB.pos;
    startLimit = Ra * sl.linkAstart + linkA.pos;
    endLimit = Ra * sl.linkAend + linkA.pos;
    c(rows) = ((sqrt(power(xb(1)-startLimit(1),2) + power(xb(2)-startLimit(2),2))) + (sqrt(power(xb(1)-endLimit(1),2) + power(xb(2)-endLimit(2),2)))) - (sqrt(power(startLimit(1)-endLimit(1),2) + power(startLimit(2)-endLimit(2),2)));
end
end

%%
function particles = updateParticles(links,particles,base,t,dt,angVel)
% Transform particle position from local to world
for i = 1 : length(particles)
	particle = particles(i);
	link = links(particle.link);
	% Compute x, the world space position of the particle.
    % Particles with relation to the base link
    R = rotationMatrix(link.angle);
    Rbase = rotationMatrix(links(base).angle);
    Rbase2 = rotationMatrix(-links(base).angle);
	x = R * particle.point + link.pos;
    xb = inv(Rbase)*x;
    Rmove = [cosd(dt*angVel) -sind(dt*angVel); sind(dt*angVel) cosd(dt*angVel)];
    Rmove2 = [cosd(0) -sind(0); sind(0) cosd(0)];
    xm = inv(rotationMatrix(t))*x;
      % INVERSE DESIGN
    %   inv(rotationMatrix(t))*x;
    %inv(Rmove)*xb;
    %Rbase2*x;
    %rotationMatrix(((length(particles)-i)/length(particles))*t*(links(base).angle))*x;
    %rotationMatrix(((length(particles)-i)/length(particles))*t*(links(base).angle))*x;
    % + (Rbase * particle.initPos - particle.initPos);
	% Append relative position to the array (grows indefinitely)
	particles(i).pointsWorld(:,end+1) = x;
    particles(i).pointsBase(:,end+1) = xb;    
    links(base).pointsA(:,end+1) = xb;    
    particles(i).pointsTry(:,end+1) = xm;
end
end

%%
function drawScene(t,links,pins,sliders,particles,angVel,pretty)

global video;

if t == 0
    clf;
    axis equal;
    hold on;
    grid on;
    xlabel('X');
    ylabel('Y');
    video = VideoWriter('output','MPEG-4');
    video.open();
end
cla;

if ~pretty
    % Draw links
    for i = 1 : length(links)
        link = links(i);
        R = rotationMatrix(link.angle);
        % Draw frame
        p = link.pos; % frame origin
        s = 0.2; % frame display size
        px = p + s*R(:,1); % frame x-axis
        py = p + s*R(:,2); % frame y-axis
        plot([p(1),px(1)],[p(2),px(2)],'r','LineWidth',3);
        plot([p(1),py(1)],[p(2),py(2)],'g','LineWidth',3);
        % Draw link geometry
        if link.grounded
            color = [0.7 0.7 0.7];
        elseif link.driver
            color = [0 1 0];
        else
            color = [0 0 1];
        end

        if link.base
            ang=0:0.01:2*pi; 
            r = link.verts(1,2)-link.verts(2,2)-0.1;
            xp=r*cos(ang);
            yp=r*sin(ang);
            plot(p(1)+xp,p(2)+yp,'LineWidth',4,'Color',[0.6,0.9,0.9]);
            length(links(i).pointsA)
            for j = 1 : length(links(i).pointsA)
                plot(links(i).pointsA(1,:),links(i).pointsA(2,:),'m');
                plot(links(i).pointsA(1,end),links(i).pointsA(2,end),'gx');
            end
        end
        E = [R,link.pos;0,0,1]; % transformation matrix
        vertsLocal = [link.verts;ones(1,size(link.verts,2))];
        vertsWorld = E * vertsLocal;
        plot(vertsWorld(1,[1:end,1]),vertsWorld(2,[1:end,1]),'Color',color);
    end
    % Draw pins
    for i = 1 : length(pins)
        pin = pins(i);
        linkA = links(pin.linkA);
        linkB = links(pin.linkB);
        Ra = rotationMatrix(linkA.angle);
        Rb = rotationMatrix(linkB.angle);
        xa = Ra * pin.pointA + linkA.pos;
        xb = Rb * pin.pointB + linkB.pos;
        plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','c');
        plot(xb(1),xb(2),'mx','MarkerSize',10,'LineWidth',2);
    end
    % Draw Sliders
    for i = 1 : length(sliders)
        sl = sliders(i);
        linkA = links(sl.linkA);
        linkB = links(sl.linkB);
        Ra = rotationMatrix(linkA.angle);
        Rb = rotationMatrix(linkB.angle);
        xa = Ra * sl.pointA + linkA.pos;
        xb = Rb * sl.pointB + linkB.pos;
        spx = Ra * sl.linkAstart + linkA.pos;
        spy = Ra * sl.linkAend + linkA.pos;
        sx = [spx(1) spy(1)]';
        sy = [spx(2) spy(2)]';
        plot(sx,sy,'LineWidth',1.4,'Color','m');
    %    plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','g');
        plot(xb(1),xb(2),'rx','MarkerSize',10,'LineWidth',2);
    end
    % Draw particles
    for i = 1 : length(particles)
        particle = particles(i);
        if ~isempty(particle.pointsWorld)
            %plot(particle.pointsWorld(1,:),particle.pointsWorld(2,:),'k');
            %plot(particle.pointsWorld(1,end),particle.pointsWorld(2,end),'ko');
            plot(particle.pointsBase(1,:),particle.pointsBase(2,:),'k');
            plot(particle.pointsBase(1,end),particle.pointsBase(2,end),'ko');
            %h = plot(particle.pointsTry(1,:),particle.pointsTry(2,:),'r');
            %hk = plot(particle.pointsTry(1,end),particle.pointsTry(2,end),'ro');
          %  rotate(h,[1 1 0],t*angVel);
        end
    end
else
    cla;
    % Draw particles
    for i = 1 : length(particles)
        particle = particles(i);
        if ~isempty(particle.pointsWorld)
            %plot(particle.pointsWorld(1,:),particle.pointsWorld(2,:),'k');
            %plot(particle.pointsWorld(1,end),particle.pointsWorld(2,end),'ko');
            plot(particle.pointsBase(1,:),particle.pointsBase(2,:),'k');
            plot(particle.pointsBase(1,end),particle.pointsBase(2,end),'ko');
            %h = plot(particle.pointsTry(1,:),particle.pointsTry(2,:),'r');
            %hk = plot(particle.pointsTry(1,end),particle.pointsTry(2,end),'ro');
          %  rotate(h,[1 1 0],t*angVel);
        end
    end
end
%axis equal;
title(sprintf('t=%.3f',t));
drawnow;

frame = getframe(gcf);
video.writeVideo(frame);
 end
