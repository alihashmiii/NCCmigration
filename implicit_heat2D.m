%% Louise Dyson 18/1/11
%% solves the diffusion equation dC/dt = ahat(d^2C/dx^2 + d^2C/dy^2)
%% on a rectangular domain given by
%% (max_x_grid,max_y_grid), the number of grid points in each direction, and (dx,dy), the spacing in each direction

% 16.05.2015 cleaned up some lines -- LJS

function C = implicit_heat2D(C,ahat,dx,dy,dt,t_steps)

[max_x_grid,max_y_grid] = size(C);

x_lower = 1;
y_lower = 1;
x_upper = max_x_grid;
y_upper = max_y_grid;

Matrix = zeros(max_x_grid*max_y_grid,max_x_grid*max_y_grid);
Number = zeros(max_x_grid,max_y_grid);
for i=1:max_x_grid
    for j=1:max_y_grid
        Number(i,j) = (j-1)*max_x_grid + i;
    end
end

%% Main equations
for i=x_lower+1:x_upper-1
    for j=y_lower+1:y_upper-1
        Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
        Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
        Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
        Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
        Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(2/dx^2 + 2/dy^2);
    end
end

%% Boundary conditions
% for j=y_lower:y_upper
%     Matrix(Number(x_lower,j),Number(x_lower,j)) = 1;
%
%     Matrix(Number(x_upper,j),Number(x_upper,j)) = 1;
% end
% for i=x_lower:x_upper
%     Matrix(Number(i,y_lower),Number(i,y_lower)) = 1;
%
%     Matrix(Number(i,y_upper),Number(i,y_upper)) = 1;
% end

% no flux BCs
for j=y_lower+1:y_upper-1
    i = x_lower;
    Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 2/dy^2);
    i = x_upper;
    Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 2/dy^2);
end
for i=x_lower+1:x_upper-1
    j = y_lower;
    Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(2/dx^2 + 1/dy^2);
    j = y_upper;
    Matrix(Number(i,y_upper),Number(i,y_upper)) = ahat/dy;
    Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
    Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
    Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(2/dx^2 + 1/dy^2);
end
i = x_lower; j = y_lower;
Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 1/dy^2);
i = x_lower; j = y_upper;
Matrix(Number(i,j),Number(i+1,j)) = - ahat*dt/dx^2;
Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 1/dy^2);
i = x_upper; j = y_lower;
Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
Matrix(Number(i,j),Number(i,j+1)) = - ahat*dt/dy^2;
Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 1/dy^2);
i = x_upper; j = y_upper;
Matrix(Number(i,j),Number(i-1,j)) = - ahat*dt/dx^2;
Matrix(Number(i,j),Number(i,j-1)) = - ahat*dt/dy^2;
Matrix(Number(i,j),Number(i,j)) = 1 + ahat*dt*(1/dx^2 + 1/dy^2);
        
Matrix = sparse(Matrix);

temp_C= reshape(C(:,:,1),max_x_grid*max_y_grid,1);
for k=1:t_steps
    temp_C = Matrix\temp_C;
    C  = reshape(temp_C,max_x_grid,max_y_grid);
end
