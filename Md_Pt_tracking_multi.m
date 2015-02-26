clear

addpath(genpath('C:\Kezhi\MyCode\.'));
load MidPt_traj

points = Frenet_Pt{1}.xy;

max_linking_distance = 30;
max_gap_closing = Inf;
debug = true;

[ tracks adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

% Dimensionality of the simulated problem (2 for 2D, 3 for 3D)
n_dim = 2;

% Number of rames to track the points over
n_frames = 57;

% Aproximative number of points per frame
n_points_per_frame = 1;

n_tracks = numel(tracks);
colors = hsv(n_tracks);

all_points = vertcat(points{:});

for i_track = 1 : n_tracks

    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.

    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);

    plot(track_points(:,2), mm-track_points(:, 1), 'Color', colors(i_track, :))

end

