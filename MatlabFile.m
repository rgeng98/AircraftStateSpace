clc; clear all; close all; format compact;
% -------------------------------------------------------------------------
%
% Bob Geng - Analyzing the longitudinal dynamics of a fighter aircraft
%
% -------------------------------------------------------------------------

%% Define the State Matrices and Calculate Transfer Function
A = [-0.5614 1 0.0717;-1.4168 -0.4932 1.6450;0 0 -10];
B = [0; 0; 10];
C = [1 0 0];
D = [0];
[num, den] = ss2tf(A,B,C,D);
H = tf(num,den)

%% Use Matlabs State Space function

state = ss(A,B,C,D);

%% Find the Controllability matrix

ctrb(state)

%% Assess the systems Controllability

if (rank(ctrb(state))==size(A,1))
    disp('Controllable');
else
    disp('Not controllable.');
end

%% Find the systems observability matrix

obsv(state)

%% Assess the systems observability

if (rank(obsv(state))==size(A,1))
    disp('Observable');
else
    disp('Not observable');
end

%% Show systems Open Loop Response

figure();
step(state)
title('Open Loop with Zero Initial Condition Step Response')

%% Show systems 0 input response

t = linspace(0, 25, 100);
r = zeros(size(t));
XrO = [0.1;0;0];
[Yr, t, Xr] = initial(state, XrO, 10);
figure();
plot(t, Yr)
title('Zero Input Response - X_0 from Step 2')

%% Compute Controller Canonical Form and Observer Canonical Form

% CCF and OCF shown at the bottom under the CCF function
[CCF, OCF] = ccf(A,B,C,D)

%% Determine the eigenvalues of the system

eig=pole(state)

%% Assess the stability of the system

stab = 0;
for i=[1:size(eig,1)]
    if real(eig(i)) < 0
        stab = stab+1;
    end
end
mstab = 0;
for i=[1:size(eig,1)]
    if real(eig(i)) == 0
        mstab = mstab+1;
    end
end

if stab == 3
    disp('Open Loop System is Asymptotically Stable')
elseif mstab ~= 0
    disp('Open Loop System is Marginably Stable')
else
    disp('Open Loop System is Unstable')
end

%% Plot Phase Diagram

figure()
plot(Xr(:,1),Xr(:,2))
xlabel('\alpha [Rad]')
ylabel('q [rad/s]')
title('Open Loop Phase Diagram')

%% Pole Placement Controller

% 1 Second Settling time
Ts=2;
Re_S=4/Ts;
% Damping of 0.7 to reduce oscillations
damping=0.7;
w=4/(damping*Ts);
i=sqrt(-1);
% Calculate desired eigenvalues
s_1=-Re_S+i*w*sqrt(1-damping^2);
s_2=-Re_S-i*w*sqrt(1-damping^2);
eigen1=s_1;
eigen2=s_2;
% define third eigenvalue to be -10 to maintain stability
eigen3=-10;
p=[eigen1 eigen2 eigen3];
% Determine K gain
K=place(A,B,p)
% Verify K gain
Kack = acker(A,B,p)
% Calculate L gain
L = place(A',C',p)'
% Verify L gain
Lack = acker(A', C', p)'
%% Test Pole Placement Control

% Create scale factor for State Feedback control
scale = -1/(C*inv(A-B*K)*B);

% Create Scale factor for Open Loop System
scaleol = -1/(C*inv(A)*B);

state1=ss(A-B*K,B,C,D);

% Plot Both Open and Closed loop responses
figure()
step(scaleol*state, scale*state1)
legend('Open Loop Response','Closed Loop Response')
title('Pole Placement Control - Normalized Step response')

%% State Observer

% Design a State Observer using the same eigenvalues as the pole placement
% control
Ar = [A -B*K;L*C A-B*K-L*C];
Br = [B;B];
Cr = [C zeros(size(C))];
Dr = D;
JbkRr = ss(Ar, Br, Cr, Dr);
t = linspace(0, 2.5, 100);
r = zeros(size(t));

% Use initial conditions provided earlier in Final
XrO = [0.1;0;0;0;0;0];
[Yr, t, Xr] = lsim(JbkRr, r, t, XrO);

% Plot Observer Error
figure();
plot(t, Xr(:,1), 'k-', t, Xr(:,4), 'b--')
legend('State', 'Observer')
title('Zero Input Observer Error')

%% Observer Compensator

Asm = [A zeros(3,1); -C 0];
Bsm = [B; 0];
Csm = [C 0];
Br = [zeros(3,1); 1];
smeig = [eigen1, eigen2, 10*real(eigen1), 15*real(eigen2)];

% Calculate a new gain matrix based on desired eigenvalues
Ksm = place(Asm, Bsm, smeig);
K = Ksm(1:3); % K being the first n places in Ke
k_b = -Ksm(4); % -ki is the last term of Ke, and so we flip the sign here

Anew = [A B*k_b -B*K; -C 0 zeros(size(C)); L*C B*k_b A-B*K-L*C];
Bnew = [Br; zeros(size(B))];
Cnew = [Csm zeros(size(C))];
sysSM = ss(Anew,Bnew,Cnew,D);

% Plot the step response of the observer based compensator
figure()
step(sysSM)
title('Observer Based Compensator Step Response')

%% LQR Controller

% Design a LQR controller
Q = [20 1 1;
     1 20 1;
     1 1 20];
R = [10];
sys = ss(A,B,C,D);

% Determine optimal K gain
[K,S,e] = lqr(sys, Q, R);
lqr_state = ss(A-B*K, B, C, D);
t = linspace(0, 6, 100);

% Calculate scale factor
scale = -1/(C*inv(A-B*K)*B);

% Plot LQR response to step input
figure()
step(scale*lqr_state, t)
title('LQR Control - Step Response')

%% Controller Canonical Form and Observer Canonical Form function

function [sysC, sysO] = ccf(A,B,C,D)

sys = ss(A,B,C,D);

Poly = poly(A);
a1 = Poly(2); a2 = Poly(3);
P = ctrb(sys);
Pc = [a2 a1  1;
      a1  1  0;
      1  0  0];
Tc = P*Pc;

Ac = inv(Tc)*A*Tc;
Bc = inv(Tc)*B;
Cc = C*Tc;

Ao = Ac';
Bo = Cc';
Co = Bc';

sysC = ss(Ac,Bc,Cc,D);
sysO = ss(Ao,Bo,Co,D);
end
