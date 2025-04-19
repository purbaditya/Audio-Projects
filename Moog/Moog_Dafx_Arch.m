function [y_temp,y_old,xd_old] = Moog_Dafx_Arch(xin,y_temp,y_old,xd_old,k,p,N)
% Function to get filter response (per sample)
h0 = p/1.3;
h1 = 0.3*h0;
h2 = 1-p;
% xd = tanh(xin-0.0001*k*y_old(end));
xd = tanh(xin - 4*k*(y_old(end) - 0.5*xin));  % Input and feedback
for i=1:N
    if i==1
        y_temp(i) = h0*xd+h1*xd_old+h2*y_old(i);
    else
        y_temp(i) = h0*y_temp(i-1)+h1*y_old(i-1)+h2*y_old(i);
    end
end
y_old = y_temp;
xd_old = xd;

end