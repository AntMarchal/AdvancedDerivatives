function [V,r_grid,t_grid] = FiniteDiff(sigma,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T,N_t,nu)

% ----------------
% Inputs:
% ----------------
%
% sigma, k, theta_Vsck : Vasicek parameter dr=k(theta-r)dt+sigma dW
% bc_left       : boundary condition on the left border,
% bc_right      : boundary condition on the right border,
% initial_cond  : @-function describing the initial condition of the problem, initial_cond = @(r) ...
%                   r is a vector (spatial grid)    
% r_min, r_max  : real values setting the extrema of the intervals in which the solution is computed
% N_r           : number of intervals in the discretization of [0,r_max]
% T             : final time of the equation
% N_t           : number of time-steps
% nu            : parameter for the theta-method time-stepping scheme. 
%                   nu=0   --> Forward Euler
%                   nu=1   --> Backward Euler
%                   nu=0.5 --> Crank-Nicholson
%
%
% ----------------
% Outputs:
% ----------------
%
% H             : matrix containing the solution of the PDE. Each column represents the solution 
%                   on the whole spatial grid at a single time
% r_grid         : spatial grid, contaiNx the nodes on which the solution V is evaluated
% t_grid        : time grid, containing the time steps at which V is computed



% set grid and grid-size.
dr = (r_max-r_min)/N_r;
r_grid = 0:dr:r_max; 
inner_grid = r_grid(2:end-1);

% set number of time steps
dt = T/N_t;
t_grid = 0:dt:T;

% init matrix containing the solution at each time step
V = zeros(length(r_grid),length(t_grid));
V(:,1)=initial_cond(r_grid)';


% for reading convenience and for computational efficiency, we define some auxiliary vectors
r = inner_grid';
r_up = inner_grid(1:end-1)'; 
r_down = inner_grid(2:end)'; 


% we want to solve
% [ I + dt*theta* A ]* V_new = V_old - dt*(1-theta) *A* V_old  + dt*F
% we need to solve the linear system
% B H_new = f_tot
% with
% B = [ I + dt*theta* A ] 
% f_tot = H_old - deltat*(1-theta) *A *H_old + dt*F

% Matrix A
A_up   = -sigma^2/(2*dr^2) - k*(theta-r_up)/(2*dr);
A_main = sigma^2/dr^2 + r;
A_down = -sigma^2/(2*dr^2) + k*(theta-r_down)/(2*dr);

A = spdiags([NaN; A_up],1,N_r-1,N_r-1)...
  + spdiags(A_main,0,N_r-1,N_r-1)...
  + spdiags([A_down; NaN],-1,N_r-1,N_r-1);

%Matrix B = [ I + deltat*theta A]

B_up   = dt*nu*(A_up);
B_main = dt*nu*(A_main) + ones(N_r-1,1);
B_down = dt*nu*(A_down);


% time-stepping loop. 

for tn = 1:N_t
    
    % Vector F
    F = [zeros(N_r-2,1);(sigma^2/(2*dr^2)...
      + k*(theta-r(end))/(2*dr)) * bc_right(tn*dt)];

    % we need to solve the linear system
    % B H_new = f_tot, with
    % B = [ I + deltat*theta A ], 
    % f_tot = H_old - deltat*(1-theta) A H_old ]

    f_tot = V(2:end-1,tn) - dt*(1-nu)*A*V(2:end-1,tn)+dt*F;
        
         
    % solve the system and store the solution
    % B is tri-diagonal, hence the system can be efficiently solved using
    % Thomas algorithm
    
    sol = thomas(B_down,B_main,B_up,f_tot);
    V(:,tn+1)=[bc_left(tn*dt);sol; bc_right(tn*dt)];

end
end
