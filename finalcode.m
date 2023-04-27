% Parameters

L = 150; % length of the square domain

t_final = 668; % simulation time

n = 50; % number of grid points in each direction

dx = L/n; % grid spacing

dt = 1; % time step size

k = 1; % thermal conductivity

T_init = 0; % initial temperature

T_source = 1; % temperature of the heat source

source_speed = 5; % speed of the heat source


% Create a mask for the heat source

[X, Y] = meshgrid((1:n)*dx, (1:n)*dx);

source_mask = ((X - 0.1*L).^2 + (Y - 0.1*L).^2 <= (0.05*L)^2);


% Initial conditions

T = T_init*ones(n,n);

T_new = T;


% Initialize variables for heat source movement

heat_source_phase = 1;



% Main loop

t = 0;

while t < t_final


 % Update temperature at interior points

 for i = 2:n-1

 for j = 2:n-1

 T_new(i,j) = T(i,j) + k*dt/(dx^2)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) - 4*T(i,j));

 end

 end



 % Apply Neumann boundary conditions at edges

 T_new(1,:) = T_new(2,:);
 T_new(n,:) = T_new(n-1,:);

 T_new(:,1) = T_new(:,2);

 T_new(:,n) = T_new(:,n-1);



 % Update temperature at heat source

 T_new(source_mask) = T_source;



 % Update heat source position
 shift = round(source_speed*dt/dx);
 if heat_source_phase == 1
 source_mask = circshift(source_mask, [0,shift]);
 if sum(sum(source_mask(:,end-shift:end))) > 0
 heat_source_phase = 2;
 end

 elseif heat_source_phase == 2
 source_mask = circshift(source_mask, [shift,-2 ]);
 if sum(sum(source_mask(end-shift:end,:))) > 0
 heat_source_phase = 3;
 end

 elseif heat_source_phase == 3
 source_mask = circshift(source_mask, [0, shift]);
 if sum(sum(source_mask(:,1:shift))) > 0
 heat_source_phase = 4;
 end

 elseif heat_source_phase == 4
 source_mask = circshift(source_mask, [0,-shift]);
 if sum(sum(source_mask(:,1:shift))) > 0
 heat_source_phase = 5;
 end
 
 elseif heat_source_phase == 5
 source_mask = circshift(source_mask, [-2,shift]);
 if sum(sum(source_mask(:,end-shift:end))) > 0
 heat_source_phase = 6;
 end
 
 else  heat_source_phase == 6;
 source_mask = circshift(source_mask, [0,-shift]);
 if sum(sum(source_mask(:,1:shift))) > 0
 heat_source_phase = 1;
 end
 end


 % Update time and temperature

 t = t + dt;

 T = T_new;



 % Plot temperature

 imagesc(T);

 colorbar;

 axis equal tight;

 title(sprintf('Time: %0.2f', t));

drawnow;

end