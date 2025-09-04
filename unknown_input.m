%  Quadratic Sum Square Error for identification of  k1,k2,...kn, alfa,beta and unknown inputs of  a n-DOF shear-building 
clear all;clc;
% time step
dt=0.01;
n=6; % The number of DDOF
l=5; % The number of acceleromete

% Measured acceleration responses
load acc.mat
accn=acc';
[ll,nn]=size(acc);

m=60*ones(1,6);                 
mass=diag(m); 
% The number of accelerometers
%dd=eye(n); % Location of Sensors
dd=zeros(l,n);
dd(1,1)=1;dd(2,2)=1;dd(3,3)=1,dd(4,4)=1;dd(5,5)=1;
% dd(6,6)=1;dd(7,7)=1;
y=dd*accn; % Measured acceleration responses

B=eye(n);
% B_un    - excitation influence matrix associated with the r-unknown excitation
Bl=[B(:,5)];      
B_un=inv(mass)*Bl;
G_un=dd*B_un;


% Initial values
X(1:2*n,1)=zeros(2*n,1);                             % Initial values of displacements and velocities
X(2*n+1:2*n+4,1)=110*ones(4,1);                      % Initial values of stiffness 
X(2*n+5,1)=1.0;                                       % Initial values of damping 
X(2*n+6,1)=0.09;

pk=zeros(2*n+3+2);
pk(1:2*n,1:2*n)=eye(2*n);                      % Initial values for error covariance of matrix 
pk(2*n+1:2*n+4,2*n+1:2*n+4)=10^8*eye(4);% Initial values of stiffness error covariance matrix
pk(2*n+5,2*n+5)=10;        % Initial values of damping error covariance matrix 
pk(2*n+6,2*n+6)=1;

Q=10^-8;                                           
% measurement noise 
R=0.001*eye(l);                                     % Recursive Solution

%load force1.txt

%f_un=force1(:,2);
 f_un(1)=0;

for k=1:ll-1;


  A=zeros(2*n+4+2);  % State transition matrix
  A(1:n,n+1:2*n)=eye(n);
  [stiff,damp,fkp,fap,fbp]=kcm(n,X(:,k));
  A(n+1:2*n,:)=-inv(mass)*[stiff,damp,fkp,fap,fbp];
  Fi=eye(2*n+4+2)+A*dt;
   
  
  hk=dd*inv(mass)*(-stiff*X(1:n,k)-damp*X(n+1:2*n,k));
  E=dd*A(n+1:2*n,:);
  
% The predicted extended state vector by numerical integration

  OPTIONS = [];
  pre=ode45(@predict,[dt*k dt*(k+1)],X(:,k),OPTIONS,n,f_un(k),B_un,mass); % Assume f_un is constant in [k*dt (k+1)*dt],dt must be small
  Xbk_1=pre.y(:,end);

  [X(:,k+1),pk_1]=klm(Xbk_1,pk,y(:,k),hk,f_un(k),Fi,E,G_un,R,Q);
   pk=pk_1;
   
   k
 [ X(2*n+1:2*n+6,k+1)]

% Estiamte the unkonw input by recusrive least squear estimation

  [stiff,damp,fkp,fap,fbp]=kcm(n,X(:,k+1));
  
  hk_1=dd*inv(mass)*(-stiff*X(1:n,k+1)-damp*X(n+1:2*n,k+1));
  [f_un(k+1)]=rlse(y(:,k+1),hk_1,G_un,R);
  
end

