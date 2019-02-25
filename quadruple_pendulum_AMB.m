% quadruple pendulum using angular momentum balance

function [t , q] = quadruple_pendulum_AMB(p,q0,reltol,abstol,tmax,timesteps)
    
    options = odeset('reltol',reltol,'abstol',abstol);
    EOMs_derive_quadruple_pendulum_AMB;
    [t , q] = ode45('quadruple_pendulum_AMB_rhs',linspace(0,tmax,timesteps),q0,options,p);

end