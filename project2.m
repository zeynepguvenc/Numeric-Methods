clear all, close all, clc
disp('Couple Dynamics:')
tol    = 0.0000001;
max_it = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.85; %cyan
b = 0.40; %red
c = 0.87; %blue
d = 0.55; %yellow
e = 0.45; %blue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetadot = 2*pi;           %% angular speed (rad/s)
period   = 2*pi/thetadot;  %% time to take for a cycle of driving arm
n  = 360;                  %% number of intervals in a cycle
cycle = 2;                 %% number of cycles to animate
N = cycle*n+1;
t  = linspace(0,cycle*period,N);  % t = time
Ax = 0;
Ay = 0;
Dx = a;
Dy = 0;
theta0 = pi/2;                       
THETA  = thetadot*t + theta0;
THETA  = mod(THETA,2*pi);    
ALPHA  = zeros(1,N);
BETA   = zeros(1,N);

triangle_angel=pi/3;
angles = [0.1, 1]';               % Initial guess for alpha and beta


for i = 1:1:N
    Bx = b*cos(THETA(i));
    By = b*sin(THETA(i));
    iter       = 1;
    while (iter <= max_it)
        alpha = angles(1);
        beta  = angles(2);

        f = [ (c*cos(alpha)-d*cos(beta)+Bx-Dx)
              (c*sin(alpha)-d*sin(beta)+By   ) ];

        J = [ (-c*sin(alpha))  ( d*sin(beta))
              ( c*cos(alpha))  (-d*cos(beta)) ]; 
        
        angles_new = angles -J\f;
        err        = norm(angles_new - angles);
        if err <= tol
            angles_next = angles_new;
            break
        else
            angles = angles_new;
        end
        iter = iter + 1;
    end
    if (iter >  max_it), error('Newton method did not converge'); end

    angles   = mod(angles_next,2*pi);    
    ALPHA(i) = angles(1);
    BETA(i)  = angles(2);
end

AX = zeros(size(t));
AY = zeros(size(t));
BX = b*cos(THETA);
BY = b*sin(THETA);
DX = a*ones(size(t));
DY = zeros(size(t));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EX = b*cos(THETA) + e*cos(triangle_angel+ALPHA);
EY = b*sin(THETA) + e*sin(triangle_angel+ALPHA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX1 = BX + c*cos(ALPHA); 
CY1 = BY + c*sin(ALPHA); 
CX2 = DX + d*cos(BETA); 
CY2 = DY + d*sin(BETA);
%%%%%%%%%%%%%%%%%%%%%%%
KX = BX + c*cos(ALPHA);
KY = BY + c*sin(ALPHA); 
MX = AX + b*cos(THETA);
MY = AY + b*sin(THETA);
%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
% ratio  = 228/127;
% factor = 0.70;
% pause(3)
for i=1:3:N
    text_plot = ['time = ',num2str(t(i),'%6.2f'),' s']; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([AX(i) BX(i)],[AY(i) BY(i)],'r-', ...
        [BX(i) CX1(i)],[BY(i) CY2(i)],'b-', ...
        [DX(i) CX2(i)],[DY(i) CY2(i)],'y-', ...
        [AX(i) DX(i)], [AY(i) DY(i)],'c-',...
        [EX(i) CX2(i)],[EY(i) CY2(i)],'b-', ...
        [BX(i) EX(i)],[BY(i) EY(i)],'b-','linewidth',3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grid on, xlabel('x'), ylabel('y'), title(text_plot), axis([-0.6 1.2 -0.6 1.2])
    %%%%%%%%%%%%%%%%%
    hold on
    plot(KX,KY,'y-')
    plot(MX,MY,'r-')
    plot(EX,EY,'b-')
    pause(0.02)
    hold off
    %%%%%%%%%%%%%%%%%
    pause(0.02)
end
figure(2)
subplot(131), plot(t,THETA,'bo','linewidth',2), grid on, xlabel('Time [s]'), ylabel('\theta [rad)')
subplot(132), plot(t,ALPHA,'ro','linewidth',2), grid on, xlabel('Time [s]'), ylabel('\alpha [rad)')
subplot(133), plot(t,BETA, 'ko','linewidth',2), grid on, xlabel('Time [s]'), ylabel('\beta [rad)')