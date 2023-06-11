setup()

%lander([3.5 1 0],[0 0.9335 0],5);

lander([3.5+0.33 1 0],[0 0.9335 0],5);

%lander([0.5 0.5 0],[0.5 0.316505 0],5);

%%

% Round 1 

setup()

launch_pos = [0.5 0.5 0]; % Initialize as launch position (somewhere on earth) 
launch_vel = 0.93; % Initialize as launch velocity
launch_angle = 29:0.2:31; % Launch angle in degrees

times = [];

for i = 1:length(launch_angle)
    vel = [launch_vel*cos(launch_angle(i)*pi/180) launch_vel*sin(launch_angle(i)*pi/180) 0];
    times(i) = lander(launch_pos,vel,2);

end

subplot(2,1,2)
plot(launch_angle,times), xlabel('angle (degrees)'), ylabel('time (s)')

%% 

% Round 3 

setup() 

global EARTH

orbiter_vel = 2;
lander_vel = 2.57;
distances = 2:0.1:3;

times = [];

for i = 1:length(distances)
    times(i) = lander([EARTH(1)+distances(i) EARTH(2) 0],[orbiter_vel lander_vel 0],2);

end

subplot(2,1,2)
plot(distances,times), xlabel('distance (m)'), ylabel('time (s)')

%%

function setup()

clear; clc; close; 

% [0 0 0] is the bottom left corner of the board

global TITAN TITAN_RAD_OUT TITAN_RAD_IN EARTH g

TITAN = [3.5 1.5];
TITAN_RAD_OUT = 0.45;  
TITAN_RAD_IN = 0.15;
EARTH = [0 0]; 

g = 9.81;

subplot(2,1,1)

% Graph gravity well 
r = linspace(0.15,0.45,30);
a = linspace(0,2*pi,30);
[R,A] = ndgrid(r,a);
Z = -(R.^2-0.913*R+0.209)./(4.93*R.^2+R+0.36);
[X,Y,Z] = pol2cart(A,R,Z);
mesh(X+TITAN(1),Y+TITAN(2),Z), hold on, grid on, xlim([0 7]), ylim([0 2.5]), zlim([-0.2 0]), axis('equal')

r = linspace(0.1272,0.15,5);
[R,A] = ndgrid(r,a);
Z = -0.457+2.02*R;
[X,Y,Z] = pol2cart(A,R,Z);
mesh(X+TITAN(1),Y+TITAN(2),Z) 

% Graph earth 
%plot3(EARTH(1),EARTH(2),0,'b.','MarkerSize',50)

end

function t = lander(pos,vel,dot_size)

global TITAN TITAN_RAD_OUT TITAN_RAD_IN g

step = 0.001; % Time step 
t = 0;

R = 0.02; % Lander radius
m = 4/3*pi*R^3*150; % Lander mass
 
rho = 1.29;
C_d = 0.47;

drag = @(vel) 1/2*rho*sqrt(vel(1)^2+vel(2)^2+vel(3)^2)*C_d*pi*R^2*vel;

while 1  

    plot3(pos(1),pos(2),pos(3),'r.','MarkerSize',dot_size)

    % Distance between lander and titan 
    r = sqrt((pos(1)-TITAN(1))^2+(pos(2)-TITAN(2))^2);
   
    if r > TITAN_RAD_OUT % Lander not in gravitaional field 
        %accel = -drag(vel)/m;
        pos = pos+vel*step;%+accel*step^2/2;
        %vel = vel+accel*step;
    
    elseif r > 0.1272
        if r > TITAN_RAD_IN % Lander in gravitational field  
            theta = atan(-1/10*(550109*r^2-134074*r-53768)/(493*r^2+100*r+36)^2); % Slope of gravity well

        else % Lander in gravitational field (linear part)
            theta = atan(2.02);
        
        end

        %[x,~,EXITFLAG] = fsolve(@(x) sys(x,m,R,r,theta,pos,vel,drag,step),[0 0 -1 1 1 1]);
        EXITFLAG = 0;

        if EXITFLAG == 1
            accel = [x(1) x(2) x(3)];
            pos = pos+vel*step+accel*step^2/2;
            vel = vel+accel*step; 
        
        else 
            accel = g*[sin(2*theta)*(TITAN(1)-pos(1))/(2*r) sin(2*theta)*(TITAN(2)-pos(2))/(2*r) -sin(theta)^2]-drag(vel)/m;
            pos = pos+vel*step+accel*step^2/2; % f(x+h) = f(x)+f'(x)*h+f''(x)*h^2/2 (second order taylor)
            vel = vel+accel*step; % Adjust velocity 
       
            % Centripetal adjustment
            norm = [TITAN(1)-pos(1) TITAN(2)-pos(2) r*tan(pi/2-theta)];
            % need atan(-vel(3)/dot([vel(1) vel(2)],[norm(1) norm(2)]/r)) = theta
            vel = vel+norm*((-vel(3)-tan(theta)/r*(vel(1)*norm(1)+vel(2)*norm(2)))/(norm(3)+tan(theta)/r*(norm(1)^2+norm(2)^2)));
    
        end

    else
        break

    end 
    
    if r > TITAN_RAD_OUT && dot([vel(1) vel(2)],[TITAN(1)-pos(1) TITAN(2)-pos(2)]) < 0
        t = 0; 
        break 

    end 

    if pos(2) > TITAN(2)-TITAN_RAD_OUT
        t = t+step; % Adjust time 

    end

    %getframe(gcf)

end

end

% x = [ax,ay,az,alpha,F_f,N]

function F = sys(x,m,R,r,theta,pos,vel,drag,step)

global TITAN g

normal = [TITAN(1)-pos(1) TITAN(2)-pos(2) r*tan(pi/2-theta)];
planar_accel = @(ax,ay,az) [ax ay az]-dot(normal/norm(normal),[ax ay az]);

F = zeros(6,1);
F(1:3) = [0 0 -m*g]-drag(vel)+x(6)*normal/norm(normal)...
    -x(5)*planar_accel(x(1),x(2),x(3))/norm(planar_accel(x(1),x(2),x(3)))-m*[x(1) x(2) x(3)];
F(4) = R*x(5)-2/5*m*R^2*x(4);
F(5) = x(4)*R-norm(planar_accel(x(1),x(2),x(3)));
F(6) = -(vel(3)+x(3)*step)/...
    dot([vel(1)+x(1)*step vel(2)+x(2)*step],[normal(1) normal(2)]/r)-tan(theta);


end
