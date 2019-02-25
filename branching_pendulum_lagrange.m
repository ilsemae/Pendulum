% branching pendulum using lagrange equations

function [t , q] = branching_pendulum_lagrange(p,q0,reltol,abstol,tmax,timesteps)
    
    options = odeset('reltol',reltol,'abstol',abstol);
    EOMs_derive_branching_pendulum_lagrange;
    [t , q] = ode45('branching_pendulum_lagrange_rhs',linspace(0,tmax,timesteps),q0,options,p);

end