function [xy,siii,ibeU]=supermeshstructurepolar(m,n)
mm=m*2;nn=n*2;
point=0;
for i=1:nn+1
      for j=1:mm+1
           point=point+1;
               xy(1,point)=(j-1)/mm;
               xy(2,point)=2*pi*(i-1)/nn;
        end
end

Le=0;%单元标号
for i=1:n
    for j=1:m
        Le=(i-1)*m+j;
         siii(1,Le)=(i-1)*2*(mm+1)+(2*j-1);
         siii(2,Le)=(i-1)*2*(mm+1)+2*j;
         siii(3,Le)=(i-1)*2*(mm+1)+2*j+1;
         siii(4,Le)=(i-1)*2*(mm+1)+(mm+1)+(2*j-1);
         siii(5,Le)=(i-1)*2*(mm+1)+(mm+1)+2*j;
         siii(6,Le)=(i-1)*2*(mm+1)+(mm+1)+2*j+1;
         siii(7,Le)=i*2*(mm+1)+(2*j-1);
         siii(8,Le)=i*2*(mm+1)+2*j;
         siii(9,Le)=i*2*(mm+1)+2*j+1;  
    end
end
% 
% ibeU1=1:2*m+1;
% ibeU2=(2*m+1)*2*n+1:(2*m+1)*(2*n+1); 
% ibeU3=(2*m+1)+1:2*m+1:(2*m+1)*(2*n-1)+1; 
% ibeU4=(2*m+1)*2:2*m+1:(2*m+1)*n*2; 
% ibeU=[ibeU1,ibeU2,ibeU3,ibeU4];
ibeU1=(2*m+1):2*m+1:(2*m+1)*(2*n+1);
ibeU=ibeU1;