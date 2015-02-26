function X = create_particles_hypo(Npop_particles,Frenet_init)

% Initialization function to generate first hypothese based on Frenet_init
%
% 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 24/02/2014




for ii = 1: Npop_particles;
    ske_spline = Frenet_init.xy + 0.3*randn(size(Frenet_init.xy));

    X{ii}.xy = ske_spline; 
    X{ii}.T = zeros(size(ske_spline));
    X{ii}.N = zeros(size(ske_spline));
    X{ii}.vel = zeros(1,2);
    X{ii}.omg = 0;
end

