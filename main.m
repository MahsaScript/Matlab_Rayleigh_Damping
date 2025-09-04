
clc
clear

h=2; 


kk=[130,100,130,130];
ElementsInformation=[ 0, 0, 1, 2;
                      1, 2, 3, 4;
                      3, 4, 5, 6;
                      5, 6, 0, 0];

stiff=zeros(6,6);
                 
 for i=1:4
         p=ElementsInformation(i,:);
 K=kk(i)*[12/h/h,   6/h,  -12/h/h,  6/h;
            6/h,     4,     -6/h,    2;
        -12/h/h,  -6/h,   12/h/h, -6/h;
            6/h,     2,     -6/h,    4;];
                 
      
          for j=1:4          
            for k=1:4
               if p(j)~=0 
                 if p(k)~=0
                    ii=p(j);
                    jj=p(k);
                    
                    stiff(ii,jj)=stiff(ii,jj)+K(j,k);
                  
                end
                  end
            end
           
          end  
      
      end

                 

                     
m=60*ones(1,6);                     
mass=diag(m);                     
                     
damp=1.1*mass+0.079*stiff;

n=6;

load force.txt;
t=force(:,1);
ff=force(:,2);



% State transition matrix
A=zeros(2*n);
A(1:n,n+1:2*n)=eye(n);
A(n+1:2*n,1:n)=-inv(mass)*stiff;
A(n+1:2*n,n+1:2*n)=-inv(mass)*damp;
   
% Excitation localization matrix 
locat=zeros(n,1); locat(5)=1.0;
B=zeros(2*n,1);
B(n+1:2*n)=inv(mass)*locat;

% initial condition
X0=zeros(2*n,1);   % initial conditions
[Y,X]=lsim(A,B,A,B,ff,t,X0);

% Acceleration response
acc=Y(:,n+1:2*n);



% % calculate the rms of all measurements
% ll=length(acc(:,1));
% for i=1:n
%    noise=randn(ll,1);
%    accn(:,i)=acc(:,i)+nl/100*std(acc(:,i))*noise;
% end;
% 
% acc=accn;

%fftplot(acc(:,3),0.01);
 ws=eig(stiff,mass); 
   f=sqrt(ws)/2/pi    
   
   save acc.mat acc
  