function X = create_particles_hypo(Npop_particles,xy)

% Initialization function to generate first hypothese based on Frenet_init
%
% 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 24/02/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.



for ii = 1: Npop_particles;
    ske_spline = xy + 0.3*randn(size(xy));

    X{ii}.xy = ske_spline; 
    X{ii}.T = zeros(size(ske_spline));
    X{ii}.N = zeros(size(ske_spline));
    X{ii}.vel = zeros(1,2);
    X{ii}.omg = 0;
    X{ii}.D = 0;
end

