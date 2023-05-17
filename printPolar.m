function printPolar(z)
% This function converts the complex numbers in the input matrix to polar
% coordinates and displays the resulting matrix in the desired format

% Get magnitude and angle of z
r = abs(z);
theta = angle(z);
theta = theta .* 180/pi;

% Convert z to polar coordinates in the desired format
z_polar = r + " < " + theta;

% Display the modified matrix in the desired format
disp(z_polar);
end
