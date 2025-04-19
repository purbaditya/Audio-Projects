function [y_temp,dy_temp_dp,dy_temp_dk,y_old,dy_old_dp,dy_old_dk,xd_old,dxd_old_dp,dxd_old_dk] = Moog_Dafx_Arch_forward_backward(xin,y_temp,dy_temp_dp,dy_temp_dk,y_old,dy_old_dp,dy_old_dk,xd_old,dxd_old_dp,dxd_old_dk,k,p,N)
% P. Bhattacharya
% Function to calculate the filter response per sample and provide the
% partial derivative value w.r.t. p (i.e. f_c) and k (TODO-i.e Q-factor?)

%% Forward - each stage: y(n) = h0*x(n)+h1*x(n-1)+h2*y(n-1)
h0 = p/1.3;
h1 = 0.3*h0;
h2 = 1-p;
% xd = tanh(xin-0.0001*k*y_old(end));         % Normal
xd = tanh(xin - 4*k*(y_old(end) - 0.5*xin));  % Input and final feedback with delay (Dafx code)
%xd = atan(xin - 4*k*(y_old(end) - 0.5*xin));
for i=1:N % N - Number of Moog stages
    if i==1
        y_temp(i) = h0*xd+h1*xd_old+h2*y_old(i);
    else
        y_temp(i) = h0*y_temp(i-1)+h1*y_old(i-1)+h2*y_old(i);
    end
end

%% Backward (before state update)
dh0 = 1/1.3;
dh1 = 0.3/1.3;
dh2 = -1;

dth = (1-xd.^2);
%dth = 1/(1+tan(xd_old).^2); % atan
dxd_dp = -4*k*dth*dy_old_dp(end);
dxd_dk = -4*dth.*(y_old(end)+k*dy_old_dk(end));

for i=1:N % N - Number of Moog stages
    if i==1
        dy_temp_dp(i) = dh0*xd+dh1*xd_old+dh2*y_old(i)+h0*dxd_dp+h1*dxd_old_dp+h2*dy_old_dp(i);
        dy_temp_dk(i) = h0*dxd_dk+h1*dxd_old_dk+h2*dy_old_dk(i);
        % TODO dy_temp_dk
    else
        dy_temp_dp(i) = dh0*y_temp(i-1)+dh1*y_old(i-1)+dh2*y_old(i)+h0*dy_temp_dp(i-1)+h1*dy_old_dp(i-1)+h2*dy_old_dp(i);
        dy_temp_dk(i) = h0*dy_temp_dk(i-1)+h1*dy_old_dk(i-1)+h2*dy_old_dk(i);
        % TODO dy_temp_dk
    end
end

%% Now update the states for next iteration
dxd_old_dp = dxd_dp;
dxd_old_dk = dxd_dk;
dy_old_dp  = dy_temp_dp;
dy_old_dk = dy_temp_dk;
y_old = y_temp;
xd_old = xd;

end