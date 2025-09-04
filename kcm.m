function [stiff,damp,fkp,fap,fbp]=kcm(n,xk)
ElementsInformation=[ 0, 0, 1, 2;
                      1, 2, 3, 4;
                      3, 4, 5, 6;
                      5, 6, 0, 0];
m=60*ones(1,6);                   
mass=diag(m);       
h=2;


% Stiffness Matrix

s1=xk(1:n);                    
s2=xk(n+1:2*n);                  
s3(1:4)=xk(2*n+1:2*n+4);         % Stiffness parameters
s4=xk(2*n+5);
s5=xk(2*n+6);


stiff=zeros(6);

 for i=1:4
         p=ElementsInformation(i,:);
 K=s3(i)*[12/h/h,   6/h,  -12/h/h,  6/h;
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
      
      %%%%%%%%%%%%%%%%%%%%%%%%
     
      damp=s4*mass+s5*stiff;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %∂‘k«Ûµº%%%%%%%%%%%%
  
  

 for i=1:4
     stiff1=zeros(6);
         p=ElementsInformation(i,:);
      K=[12/h/h,   6/h,  -12/h/h,  6/h;
            6/h,     4,     -6/h,    2;
        -12/h/h,  -6/h,   12/h/h, -6/h;
            6/h,     2,     -6/h,    4;];
                 
      
          for j=1:4          
            for k=1:4
               if p(j)~=0 
                 if p(k)~=0
                    ii=p(j);
                    jj=p(k);
                    stiff1(ii,jj)=K(j,k);
                  
                end
                  end
              end
          end  
                  fkp1(:,i)=stiff1*s1;
                  fkp2(:,i)=stiff1*s2;
 
      end
          
              

    fkp=fkp1+fkp2*s5;

    fap=mass*s2;
    fbp=stiff*s2;               
         
    
    
    
    
    
    