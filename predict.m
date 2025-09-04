function yp=predict(tb,yb,n,f_unk,B_un,mass);

% functions of the generalized state equations of a Linear System

g(1:n,1)=yb(n+1:2*n);
[stiff,damp,fkp,fap,fbp]=kcm(n,yb);
g(n+1:2*n)=B_un*f_unk+inv(mass)*(-stiff*yb(1:n)-damp*yb(n+1:2*n));
g(2*n+1:2*n+4+2)=zeros(4+2,1);

yp =g;

   