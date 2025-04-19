function [x,m,v,g] = optimize(x,dydx,m,v,g,n,lr_multiplier,scale,layer,f_update,updaterule,optimizeall)
if ~optimizeall && layer == f_update   % update a specific layer per iter
    if strcmp(updaterule,'sgd')==1            % stochastic gradient descent
        lr = 1*lr_multiplier;
        x = x - lr*dydx(n);
    else                             % adaptive stochastic gradient methods 
        % parameters %
        b1 = 0.9;
        b2 = 0.999;
        lr = 1*lr_multiplier;
        % momentum (m), velocity (v) and learning rate (lr) updated/adapted
        m = b1*m + (1-b1)*dydx(n);
        if strcmp(updaterule,'adagrad')==0
            v = b2*v + (1-b2)*dydx(n)^2;
        else
            v = v + dydx(n)^2;
        end
        mh = m./(1-b1^n);
        vh = v./(1-b2^n);
        if strcmp(updaterule,'diffgrad')==1
            gd = 1./(1+exp(-abs(g-dydx(n))));
            x = x - gd.*lr*(mh/(sqrt(vh)+eps));
        elseif strcmp(updaterule,'adam')==1
            x = x - lr*(mh/(sqrt(vh)+eps));
        elseif strcmp(updaterule,'momentum')==1
            x = x - lr*m;
        elseif strcmp(updaterule,'rmsprop')==1
            x = x - lr*dydx(n)/(sqrt(v)+eps);
        elseif strcmp(updaterule,'adagrad')==1
%             if x < 0
%                 x = min(-0.1,x - lr*200*dydx(n)/(sqrt(v)+eps));
%             else
%                 x = max(0.1,x - lr*200*dydx(n)/(sqrt(v)+eps));
%             end
            x = x - lr*(scale*dydx(n)/(sqrt(v)+eps));
        end
        g = dydx(n);
    end
else % update all layers (filters) per iter
    if strcmp(updaterule,'sgd')==1
        lr = 1*lr_multiplier;
        x = x - lr*dydx(n);
    else
        % parameters %
        b1 = 0.9;
        b2 = 0.999;
        lr = 1*lr_multiplier;
        % momentum (m), velocity (v) and learning rate (lr) updated/adapted
        m = b1*m + (1-b1)*dydx(n);
        if strcmp(updaterule,'adagrad')==0
            v = b2*v + (1-b2)*dydx(n)^2;
        else
            v = v + dydx(n)^2;
        end
        mh = m./(1-b1^n);
        vh = v./(1-b2^n);
        if strcmp(updaterule,'diffgrad')==1
            gd = 1./(1+exp(-abs(g-dydx(n))));
            x = x - gd.*lr*(mh/(sqrt(vh)+eps));
        elseif strcmp(updaterule,'adam')==1
            x = x - lr*scale*(mh/(sqrt(vh)+eps));
        elseif strcmp(updaterule,'momentum')==1
            x = x - lr*m;
        elseif strcmp(updaterule,'rmsprop')==1
            x = x - lr*dydx(n)/(sqrt(v)+eps);
        elseif strcmp(updaterule,'adagrad')==1
            x = x - lr*(scale*dydx(n)/(sqrt(v)+eps));
        end
        g = dydx(n);
    end
end
