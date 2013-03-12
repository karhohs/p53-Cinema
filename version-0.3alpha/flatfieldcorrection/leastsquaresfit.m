%----- SUBFUNCTION LEASTSQUARESFIT -----
function [a,b]=leastsquaresfit(x,y,j,k)
xm=mean(x);
ym=mean(y);
SSxx=sum(x.*x)-length(x)*xm^2;
SSyy=sum(y.*y)-length(y)*ym^2;
SSxy=sum(x.*y)-length(x)*xm*ym;
b=SSxy/SSxx;
a=ym-b*xm;
r2=(SSxy^2)/(SSxx*SSyy);
if r2<0.8
    disp(['not a good fit. x=', num2str(k), ' y=', num2str(j)])
end
end