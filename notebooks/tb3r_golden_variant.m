clearvars; clc;

% Delete any P-code file for this script to ensure changes take effect
pcode_file = [mfilename('fullpath'), '.p'];
if exist(pcode_file, 'file')
    delete(pcode_file);
    fprintf('Deleted obsolete P-code file: %s\n', pcode_file);
end

%% Setup OceanMesh2D paths
run('/home/pj24001722/ku40000343/Github/OceanMesh2D/setup_oceanmesh2d.m');
addpath('/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/data');  % Add Tokyo Bay regional data to path
addpath('/home/pj24001722/ku40000343/Github/OceanMesh2D/utilities/shapefile');  % Add shapefile utilities for safe_geodata
addpath('/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/data/Futtsu_coastline');  % Add Futtsu coastline subdirectory


%% MESH GENERATION CONFIGURATION
% Maximum iterations for mesh generation (affects quality vs computation time)
% Higher values = better quality but longer computation
% Typical values: 50 (fast), 100 (balanced), 200 (high quality)
MESH_ITMAX = 50;  % Use 100 for balanced quality and performance

% Validate itmax value
if MESH_ITMAX < 10 || MESH_ITMAX > 500
    warning('MESH_ITMAX=%d is outside typical range [10-500]. Mesh quality may be affected.', MESH_ITMAX);
end

fprintf('Mesh generation will use itmax=%d iterations\n', MESH_ITMAX);

%% Plotting control flag
% Set this to control figure generation and saving
% true  = Generate and save all figures (default)
% false = Skip all plotting and figure saving (for batch processing or headless servers)
ENABLE_PLOTTING = false;  % Change to false to disable all plotting

%% Configure display settings based on plotting flag
if ENABLE_PLOTTING
    fprintf('Plotting ENABLED - figures will be generated and saved\n');
    plot_setting = 1;  % Enable interactive plotting during mesh generation
else
    fprintf('Plotting DISABLED - skipping all figure generation\n');
    plot_setting = 0;  % Disable interactive plotting during mesh generation
    set(groot, 'defaultFigureVisible', 'off');
    try
        opengl software;  % Use software rendering
    catch
        % OpenGL might not be available in true headless mode
    end
end

% Create output directories if they don't exist
if ENABLE_PLOTTING && ~exist('../outputs/PNG', 'dir')
    mkdir('../outputs/PNG');
end
if ~exist('/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/outputs/meshes', 'dir')
    mkdir('/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/outputs/meshes');
end

% Setup output directories
outdir = '/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/outputs/PNG';
meshdir = '/home/pj24001722/ku40000343/Github/OceanMesh2D/Tokyo_Bay/outputs/meshes';


%% STEP 1: set mesh extents and set parameters for mesh.
x_bond_01 = [157.4634361 158.84262217 160.046858 161.07987092 ...
    161.94757837 162.65742188 163.21770287 163.63704393 ...
    163.92386584 164.08625707 164.1323057 164.06882521 ...
    163.9009549 163.63627547 163.31281843 162.89825466 ...
    162.38857843 161.80134206 161.13903484 160.4052709 ...
    159.60336662 158.7360147 157.80598566 156.81590455 ...
    155.7687459 154.66746366 153.5151905 152.3154099 ...
    151.0718658 149.78866673 148.4702309 147.12135316 ...
    145.74715292 144.35301788 142.94462003 141.52778842 ...
    140.10846886 138.69263968 137.2862313 135.89504337 ...
    134.52470652 133.18055169 131.86766327 130.59072034 ...
    129.35412832 128.16190388 127.01767503 125.92498131 ...
    124.88679107 123.90615375 122.98581257 122.12855282 ...
    121.60387397 115.64257669 158.59463098 167.27851682 ...
    157.4634361];
y_bond_01 = [51.36354136 50.35540912 49.27890879 48.14482943 ...
    46.96329793 45.74371019 44.4946652 43.2239559 ...
    41.93859546 40.64492778 39.34870289 38.0551143 ...
    36.76910533 35.56167205 34.36768078 33.12965026 ...
    31.91772707 30.7320848 29.57667543 28.45502315 ...
    27.37062929 26.32729387 25.32881846 24.37921688 ...
    23.48234609 22.64223232 21.86295783 21.14846084 ...
    20.50262818 19.92912188 19.4314301 19.01266254 ...
    18.67552551 18.42237308 18.25486888 18.17421142 ...
    18.18096748 18.27508101 18.45587251 18.7220402 ...
    19.07182143 19.5027674 20.01215433 20.59669911 ...
    21.25297733 21.97724182 22.76552456 23.61401919 ...
    24.51851667 25.47512814 26.47984282 27.52881366 ...
    28.22438211 39.405273 75.05163401 62.25276725 ...
    51.36354136];
bbox_01      = [x_bond_01',y_bond_01'];
min_el_01    = 10e3;  	   % minimum resolution in meters.
max_el_01    = 50e3; 	   % maximum resolution in meters.
wl_01        = 30;         % 60 elements resolve M2 wavelength.
dt_01        = 0;          % Automatically set timestep based on nearshore res
grade_01     = 0.3;       % mesh grade in decimal percent.
R_01         = 3; 		   % Number of elements to resolve feature.
slp_01       = 50;         % 2*pi/number of elements to resolve slope
fl_01        = -50;        % use filter equal to Rossby radius divided by 50
dem_srtm = '/home/pj24001722/ku40000343/Github/OceanMesh2D/datasets/TokyoBay/dem/SRTM15_pacific_4min.nc'; dem_kanto = '/home/pj24001722/ku40000343/Github/OceanMesh2D/datasets/TokyoBay/dem/SRTM15_kanto_15s.nc';
coastline_gshhs = 'GSHHS_f_L1';
gdat_01 = safe_geodata('shp',coastline_gshhs,...
    'dem',dem_srtm,...
    'h0',min_el_01,...
    'bbox',bbox_01);
fh_01 = edgefx('geodata',gdat_01,...
    'fs',R_01,'wl',wl_01,...
    'slp',slp_01,'fl',fl_01,...
    'max_el',max_el_01,...
    'dt',dt_01,'g',grade_01);
%% STEP 2: specify geographical datasets and process the geographical data
x_bond_02 = [138.196907945361,139.784782123012,141.294666567004,...
    139.718998032166,138.285893837899,138.196907945361];
y_bond_02 = [34.6003740729489,33.5207179317026,35.2698037998702,...
    35.9357551145740,35.2841759386230,34.6003740729489];

x_bond_03 = [139.142156 139.142156	139.144724 139.152351 139.164815 ...
    139.18178	139.202839	139.227552	139.25548	139.286203 ...
    139.31934	139.354542	139.391505	139.429958	139.469663 ...
    139.510408	139.552007	139.594291	139.637105	139.680308 ...
    139.723765	139.767348	139.810932	139.854389	139.897592 ...
    139.940406	139.98269	140.024289	140.065034	140.104739 ...
    140.143192	140.180155	140.215357	140.248493	140.279217 ...
    140.307145	140.331858	140.352917	140.369882	140.382346 ...
    140.389973	140.392541	140.392541	139.142156];
y_bond_03 = [35.859592	 35.225989	35.182486  35.139579  35.097818 ...
    35.057672  35.019512  34.983607  34.95014  34.919216 ...
    34.890891  34.865177  34.842062  34.82152  34.803513 ...
    34.788004  34.774954  34.764328  34.756096 34.750234 ...
    34.746724  34.745556  34.746724  34.750234 34.756096 ...
    34.764328  34.774954  34.788004  34.803513 34.82152 ...
    34.842062  34.865177  34.890891  34.919216 34.95014	...
    34.983607  35.019512  35.057672  35.097818 35.139579 ...
    35.182486  35.225989  35.859592	 35.859592];
%% STEP 3: create an edge function class
% set bbox02
bbox_02       = [x_bond_02',y_bond_02'];
min_el_02    = 1e3;
max_el_02    = 2e3;
wl_02        = 30;
dt_02        = 0;
grade_02     = 0.1;
R_02         = 3;
slp_02       = 50;
fl_02        = -50;

coastline_02 = '/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/outputs/tb_varres_3r/land_osm_wide';
gdat_02 = safe_geodata('shp',coastline_02,...
    'dem', dem_kanto, ...
    'h0',min_el_02,...
    'bbox',bbox_02);
fh_02 = edgefx('geodata',gdat_02,...
    'fs',R_02,'wl',wl_02,...
    'slp',slp_02,'fl',fl_02,...
    'max_el',max_el_02,...
    'dt',dt_02,'g',grade_02);

% set bbox03
bbox_03       = [x_bond_03',y_bond_03'];
min_el_03    = 1e2;
max_el_03    = 5e2;
wl_03        = 30;
dt_03        = 0;
grade_03     = 0.1;
R_03         = 3;
slp_03       = 50;
fl_03        = -50;

coastline_03 = 'coastline_2';
gdat_03 = safe_geodata('shp',coastline_02,...
    'dem', dem_kanto, ...
    'h0',min_el_03,...
    'bbox',bbox_03);
fh_03 = edgefx('geodata',gdat_03,...
    'fs',R_03,'wl',wl_03,...
    'slp',slp_03,'fl',fl_03,...
    'max_el',max_el_03,...
    'dt',dt_03,'g',grade_03);
%% STEP 4: Pass your edgefx class object along with some meshing options
mshopts = meshgen(...
    'ef',{fh_01 fh_02 fh_03},...
    'bou',{gdat_01 gdat_02 gdat_03},...
    'plot_on',plot_setting,...
    'proj','trans',...
    'itmax',MESH_ITMAX);
mshopts = mshopts.build;
% Save mesh generation progress
if ENABLE_PLOTTING
    fprog = figure('Visible','off'); 
    axp = axes('Parent',fprog); 
    hold(axp,'on');
    axes(axp);
    h = triplot(mshopts.grd.t, mshopts.grd.p(:,1), mshopts.grd.p(:,2), 'k-');
    set(h, 'Color', 'k');
    axis(axp,'equal'); 
    axis(axp,'tight');
    title(axp,'Mesh generation progress (final mesh)');
    drawnow;
    mesh_prog_file = fullfile(outdir,'mesh_generation_progress.png');
    if exist('exportgraphics','file')
        exportgraphics(axp, mesh_prog_file, 'Resolution',150);
    else
        print(fprog, mesh_prog_file, '-dpng', '-r150');
    end
    close(fprog);
end

% Close mesh generation figures
if ENABLE_PLOTTING
    all_figs = get(0, 'Children');
    for fig = all_figs'
        if ishandle(fig)
            close(fig);
        end
    end
else
    % Close all figures without saving
    close all;
end 

%% Plot and save the msh class object/write to fort.14
m = mshopts.grd; % get out the msh object
m = interp(m,{gdat_01 gdat_02 gdat_03},'mindepth',0.01);
m = make_bc(m,'auto',gdat_01,'both');

%% Save mesh data
mesh_name = 'tb_varres_3regions';
save(fullfile(meshdir, mesh_name),'m');  
write(m, fullfile(meshdir, mesh_name));
fprintf('Mesh saved as: %s\n', mesh_name);

%% Generate PNG files for visualization
if ENABLE_PLOTTING
    fprintf('\nGenerating visualization plots...\n');
else
    fprintf('\nSkipping PNG generation (plotting disabled)\n');
end

%% Figure 1: Mesh boundaries with region boxes (FULL VIEW)
if ENABLE_PLOTTING
    fig1 = findobj('Type','figure','Number',1);
    if ~isempty(fig1)
        close(fig1);
    end
    f1 = figure(1); clf(f1);
set(f1,'Name','Mesh Boundaries with Region Boxes','NumberTitle','off','Visible','on');
ax1 = axes('Parent',f1); hold(ax1,'on');

% Plot region boxes FIRST and keep their handles
r1 = plot(ax1, [x_bond_01 x_bond_01(1)], [y_bond_01 y_bond_01(1)], 'b--','LineWidth',2.5);
r2 = plot(ax1, [x_bond_02 x_bond_02(1)], [y_bond_02 y_bond_02(1)], 'g--','LineWidth',2.5);
r3 = plot(ax1, [x_bond_03 x_bond_03(1)], [y_bond_03 y_bond_03(1)], 'c--','LineWidth',2.5);

% Set viewport from boxes and lock it
all_x = [x_bond_01, x_bond_02, x_bond_03];
all_y = [y_bond_01, y_bond_02, y_bond_03];
x_margin = 0.5; y_margin = 0.5;
xlim(ax1, [min(all_x)-x_margin, max(all_x)+x_margin]);
ylim(ax1, [min(all_y)-y_margin, max(all_y)+y_margin]);
axis(ax1,'equal'); axis(ax1,'manual');

% Plot mesh boundaries AFTER boxes on the same axes handle
if ~isempty(m.op)
    for nb = 1:m.op.nope
        xop = m.p(m.op.nbdv(1:m.op.nvdll(nb),nb),1);
        yop = m.p(m.op.nbdv(1:m.op.nvdll(nb),nb),2);
        plot(ax1, xop, yop, 'b-','LineWidth',1.2);
    end
end
if ~isempty(m.bd)
    for nb = 1:m.bd.nbou
        xbd = m.p(m.bd.nbvv(1:m.bd.nvell(nb),nb),1);
        ybd = m.p(m.bd.nbvv(1:m.bd.nvell(nb),nb),2);
        if m.bd.ibtype(nb)==20  % Land boundary
            plot(ax1, xbd, ybd, 'k-','LineWidth',1.2);
        else
            plot(ax1, xbd, ybd, 'g-','LineWidth',1.2);
        end
    end
end

% FORCE boxes on top using their handles
uistack([r1 r2 r3], 'top');

xlim(ax1, [min(all_x)-x_margin, max(all_x)+x_margin]);
ylim(ax1, [min(all_y)-y_margin, max(all_y)+y_margin]);
xlabel(ax1,'Longitude'); ylabel(ax1,'Latitude');
title(ax1,'Tokyo Bay 3 Regions (Variable Resolution) - Mesh Boundaries with Region Boxes');
drawnow;

    fig1_file = fullfile(outdir,'3regions_varres_mesh_boundary.png');
    print(f1, fig1_file, '-dpng', '-r300');
    fprintf('  - Saved: 3regions_varres_mesh_boundary.png (Full mesh with boundaries)\n');
end

%% Figure 2: DEM plot for main region
if ENABLE_PLOTTING
    figure(2);
clf;
% Create DEM visualization without using geodata plot method
try
    % Create a grid for DEM visualization
    stride = gdat_01.h0/111e3;  % Convert to degrees
    xx = gdat_01.bbox(1,1):stride:gdat_01.bbox(1,2);
    yy = gdat_01.bbox(2,1):stride:gdat_01.bbox(2,2);
    [demx, demy] = meshgrid(xx, yy);
    
    % Interpolate bathymetry onto grid
    demz = gdat_01.Fb(demx, demy);
    
    % Plot using pcolor
    pcolor(demx, demy, demz);
    shading flat;
    colormap(jet);
    colorbar;
    axis equal;
    xlabel('Longitude');
    ylabel('Latitude');
    title('Tokyo Bay Region - DEM');
catch
    % Fallback if DEM data not available
    text(0.5, 0.5, 'DEM data not available', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    axis([0 1 0 1]);
    axis off;
    title('Tokyo Bay Region - DEM (No data)');
end
    print(gcf, fullfile(outdir, '3regions_varres_dem_main.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_dem_main.png\n');
end

%% Figure 3: Bathymetry plots with mesh and without
% 3a. Bathymetry plot with colormap and mesh
if ENABLE_PLOTTING
    figure(3);
clf; % Clear figure
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',m.b,'FaceColor','flat','EdgeColor','black','LineWidth',0.1);
axis equal;
axis tight; % Fit axis to data
xlim([min(m.p(:,1)) max(m.p(:,1))]);
ylim([min(m.p(:,2)) max(m.p(:,2))]);
title('Tokyo Bay Variable Resolution - Bathymetry with Mesh (m)');
xlabel('Longitude');
ylabel('Latitude');
c = colorbar;
ylabel(c, 'Depth (m)');
colormap(jet);
drawnow; % Force drawing
    print(gcf, fullfile(outdir, '3regions_varres_bathymetry_mesh.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_bathymetry_mesh.png\n');
end

% 3b. Bathymetry plot with colormap only (no mesh)
if ENABLE_PLOTTING
    figure(31);
clf;
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',m.b,'FaceColor','flat','EdgeColor','none');
axis equal;
axis tight;
xlim([min(m.p(:,1)) max(m.p(:,1))]);
ylim([min(m.p(:,2)) max(m.p(:,2))]);
title('Tokyo Bay Variable Resolution - Bathymetry (m)');
xlabel('Longitude');
ylabel('Latitude');
c = colorbar;
ylabel(c, 'Depth (m)');
colormap(jet);
drawnow;
    print(gcf, fullfile(outdir, '3regions_varres_bathymetry.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_bathymetry.png\n');
end

%% Figure 4: Resolution plots
if ENABLE_PLOTTING
    fprintf('  Calculating mesh resolution...\n');
else
    fprintf('  Skipping resolution calculation (plotting disabled)\n');
end

if ENABLE_PLOTTING
% Calculate resolution at each node (in meters)
res_node = nan(size(m.p,1),1);
for i = 1:size(m.p,1)
    % Find elements containing this node
    [rows,~] = find(m.t == i);
    if ~isempty(rows)
        % Calculate average element size around this node
        sizes = zeros(length(rows),1);
        for j = 1:length(rows)
            v = m.p(m.t(rows(j),:),:);
            % Calculate edge lengths in meters (approximate)
            edge_lengths = zeros(3,1);
            for k = 1:3
                p1 = v(k,:);
                p2 = v(mod(k,3)+1,:);
                % Approximate distance in meters
                edge_lengths(k) = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2) * 111000; % rough conversion
            end
            sizes(j) = mean(edge_lengths);
        end
        res_node(i) = mean(sizes);
    end
end
% Fill any NaN values
if any(isnan(res_node))
    res_node(isnan(res_node)) = mean(res_node(~isnan(res_node)));
end

    % 4a. Resolution plot with mesh
    figure(4);
clf;  % Clear the figure first
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',res_node,'FaceColor','flat','EdgeColor','black','LineWidth',0.1);
title('Tokyo Bay Variable Resolution - Mesh Resolution with Mesh (m)');
colorbar;
colormap(jet);
axis equal;
axis tight;
xlim([min(m.p(:,1)) max(m.p(:,1))]);
ylim([min(m.p(:,2)) max(m.p(:,2))]);
xlabel('Longitude');
ylabel('Latitude');
drawnow;  % Force drawing
    print(gcf, fullfile(outdir, '3regions_varres_resolution_mesh.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_resolution_mesh.png\n');

    % 4b. Resolution plot without mesh
    figure(41);
clf;  % Clear the figure first
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',res_node,'FaceColor','flat','EdgeColor','none');
title('Tokyo Bay Variable Resolution - Mesh Resolution (m)');
colorbar;
colormap(jet);
axis equal;
axis tight;
xlim([min(m.p(:,1)) max(m.p(:,1))]);
ylim([min(m.p(:,2)) max(m.p(:,2))]);
xlabel('Longitude');
ylabel('Latitude');
drawnow;  % Force drawing
    print(gcf, fullfile(outdir, '3regions_varres_resolution.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_resolution.png\n');
end

%% Figure 5: Full mesh plot
if ENABLE_PLOTTING
    figure(5);
clf;
triplot(m.t, m.p(:,1), m.p(:,2), 'k-', 'LineWidth', 0.1);
axis equal;
title('Tokyo Bay Variable Resolution - Full Mesh');
xlabel('Longitude');
ylabel('Latitude');
    print(gcf, fullfile(outdir, '3regions_varres_full_mesh.png'), '-dpng', '-r300');
    fprintf('  - Saved: 3regions_varres_full_mesh.png\n');
end

%% Figure 6: Mesh quality analysis
if ENABLE_PLOTTING
    fprintf('  Calculating mesh quality metrics...\n');
else
    fprintf('  Skipping mesh quality analysis (plotting disabled)\n');
end

% Calculate mesh quality metrics (always needed for statistics)
points = m.p;  
triangles = m.t;
% Use the OceanMesh2D built-in triangleAngles function
angles = triangleAngles(points, triangles); 
min_angle_per_tri = min(angles, [], 2);

if ENABLE_PLOTTING

    %% Figure 7: Histogram of minimum angles
    figure(6);
histogram(min_angle_per_tri);
xlabel('Minimum angle (degrees)');
ylabel('Number of triangles');
title('Mesh Quality: Distribution of Minimum Triangle Angles');
    print(gcf, fullfile(outdir, '3regions_varres_min_angle_histogram.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_min_angle_histogram.png\n');

    %% Figure 8: Pie chart of angle distribution
    edges  = [0, 20, 30, 40, 180];
    labels = {'< 20', '20-30', '30-40', '> 40'};
    counts = histcounts(min_angle_per_tri, edges);
    
    figure(7);
pie(counts, labels);
title('Minimum Triangle Angle Distribution (Pie Chart)');
    print(gcf, fullfile(outdir, '3regions_varres_min_angle_pie.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_min_angle_pie.png\n');
end

%% Figure 9: Variable resolution demonstration
if ENABLE_PLOTTING
    figure(9);
clf;
% Focus on transition zones to show variable resolution
trans_x = [139.0, 140.5];
trans_y = [34.5, 36.0];

% Find mesh elements in this region
in_region = m.p(:,1) >= trans_x(1) & m.p(:,1) <= trans_x(2) & ...
            m.p(:,2) >= trans_y(1) & m.p(:,2) <= trans_y(2);
region_nodes = find(in_region);

if ~isempty(region_nodes)
    % Find triangles that have at least one node in the region
    region_tri_mask = any(ismember(m.t, region_nodes), 2);
    region_tri = m.t(region_tri_mask, :);
    
    triplot(region_tri, m.p(:,1), m.p(:,2), 'b-', 'LineWidth', 0.8);
    xlim(trans_x); ylim(trans_y);
    axis equal;
    title('Tokyo Bay Variable Resolution - Transition Zones Detail');
    xlabel('Longitude'); ylabel('Latitude');
else
    text(0.5, 0.5, 'No mesh data in transition region', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    axis([0 1 0 1]); axis off;
    title('Tokyo Bay Variable Resolution - Transition Zones (No data)');
end
    print(gcf, fullfile(outdir, '3regions_varres_transition_detail.png'), '-dpng', '-r300');
    fprintf('  - Saved: 3regions_varres_transition_detail.png\n');
end

%% Figure 10: Boundary condition visualization
if ENABLE_PLOTTING
    figure(10);
clf;
if ~isempty(m.bd)
    hold on;
    for nb = 1:m.bd.nbou
        xbd = m.p(m.bd.nbvv(1:m.bd.nvell(nb),nb),1);
        ybd = m.p(m.bd.nbvv(1:m.bd.nvell(nb),nb),2);
        if m.bd.ibtype(nb)==20  % Land boundary
            plot(xbd, ybd, 'k-','LineWidth',1.5);
        elseif m.bd.ibtype(nb)==0  % Open boundary
            plot(xbd, ybd, 'r-','LineWidth',2);
        else
            plot(xbd, ybd, 'g-','LineWidth',1.5);
        end
    end
    axis equal; axis tight;
    title('Tokyo Bay Variable Resolution - Boundary Conditions');
    xlabel('Longitude'); ylabel('Latitude');
    legend({'Land', 'Open', 'Other'}, 'Location', 'best');
else
    text(0.5, 0.5, 'No boundary data available', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    axis([0 1 0 1]); axis off;
    title('Tokyo Bay Variable Resolution - Boundary Conditions (No data)');
end
    print(gcf, fullfile(outdir, '3regions_varres_boundaries.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_boundaries.png\n');
end

%% Figure 11: Element size distribution
if ENABLE_PLOTTING
    figure(11);
clf;
% Calculate element areas
element_areas = zeros(size(m.t,1),1);
for i = 1:size(m.t,1)
    v = m.p(m.t(i,:),:);
    % Area using cross product formula
    element_areas(i) = 0.5 * abs((v(2,1)-v(1,1))*(v(3,2)-v(1,2)) - (v(3,1)-v(1,1))*(v(2,2)-v(1,2)));
end

histogram(log10(element_areas));
xlabel('Log10(Element Area) [degrees^2]');
ylabel('Number of elements');
title('Tokyo Bay Variable Resolution - Element Size Distribution');
    print(gcf, fullfile(outdir, '3regions_varres_element_sizes.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_element_sizes.png\n');
end

%% Figure 12: Depth contours
if ENABLE_PLOTTING
    figure(12);
clf;
try
    % Create a regular grid for contour plotting
    x_range = linspace(min(m.p(:,1)), max(m.p(:,1)), 300);
    y_range = linspace(min(m.p(:,2)), max(m.p(:,2)), 300);
    [X, Y] = meshgrid(x_range, y_range);
    
    % Interpolate bathymetry to regular grid
    F = scatteredInterpolant(m.p(:,1), m.p(:,2), m.b, 'natural', 'none');
    Z = F(X, Y);
    
    % Define depth levels for contours
    depth_levels = [-100:5:0];  % Every 5m from -100m to 0m
    
    % Create filled contour plot
    [C, h] = contourf(X, Y, Z, depth_levels);
    
    % Add contour lines on top
    hold on;
    [C2, h2] = contour(X, Y, Z, depth_levels, 'k-', 'LineWidth', 0.3);
    
    % Add labels to major contour lines
    major_depths = [-100:20:0];  % Label every 20m
    [C3, h3] = contour(X, Y, Z, major_depths, 'k-', 'LineWidth', 1);
    clabel(C3, h3, 'LabelSpacing', 300, 'FontSize', 8, 'Color', 'k');
    
    % Customize colormap and appearance
    colormap(flipud(parula));  % Flip colormap for depth (blue=deep)
    c = colorbar;
    ylabel(c, 'Depth (m)');
    caxis([min(m.b) 0]);  % Set color limits
    
    axis equal;
    axis tight;
    title('Tokyo Bay Variable Resolution - Depth Contours (m)');
    xlabel('Longitude');
    ylabel('Latitude');
    hold off;
catch ME
    % Fallback: show bathymetry as patch if contours fail
    fprintf('  Warning: Contour generation failed, using patch visualization\n');
    patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',m.b,'FaceColor','flat','EdgeColor','none');
    colormap(flipud(parula));
    c = colorbar;
    ylabel(c, 'Depth (m)');
    axis equal; axis tight;
    title('Tokyo Bay Variable Resolution - Bathymetry (m)');
    xlabel('Longitude'); ylabel('Latitude');
end
    print(gcf, fullfile(outdir, '3regions_varres_depth_contours.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_depth_contours.png\n');
end

%% Figure 13: Regional mesh statistics
if ENABLE_PLOTTING
    figure(13);
clf;
% Create bar chart of mesh statistics per region
region_names = {'Pacific', 'Mid Resolution', 'High Resolution'};
try
    % Approximate node counts per region based on resolution
    % Rough estimates based on bounding boxes
    nodes_r1 = sum(m.p(:,1) < 138.5);  % Pacific nodes
    nodes_r2 = sum(m.p(:,1) >= 138.5 & m.p(:,1) < 139.8);  % Mid resolution
    nodes_r3 = sum(m.p(:,1) >= 139.8);  % High resolution
    
    node_counts = [nodes_r1, nodes_r2, nodes_r3];
    bar(node_counts);
    set(gca, 'XTickLabel', region_names);
    ylabel('Number of nodes');
    title('Tokyo Bay Variable Resolution - Nodes per Region');
    for i = 1:length(node_counts)
        text(i, node_counts(i)+max(node_counts)*0.02, num2str(node_counts(i)), ...
             'HorizontalAlignment', 'center');
    end
catch
    text(0.5, 0.5, 'Could not calculate regional statistics', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    axis([0 1 0 1]); axis off;
    title('Tokyo Bay Variable Resolution - Regional Statistics (Error)');
end
    print(gcf, fullfile(outdir, '3regions_varres_regional_stats.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_regional_stats.png\n');
end

%% Figure 14: Overall summary plot
if ENABLE_PLOTTING
    figure(14);
clf;
% Create 2x2 subplot with key visualizations
subplot(2,2,1);
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',m.b,'FaceColor','flat','EdgeColor','none');
axis equal; axis tight;
title('Bathymetry');
colorbar;

subplot(2,2,2);
triplot(m.t, m.p(:,1), m.p(:,2), 'k-', 'LineWidth', 0.1);
axis equal; axis tight;
title('Mesh Structure');

subplot(2,2,3);
histogram(min_angle_per_tri);
xlabel('Min angle (deg)');
ylabel('Count');
title('Quality Distribution');

subplot(2,2,4);
patch('Faces',m.t,'Vertices',m.p,'FaceVertexCData',res_node,'FaceColor','flat','EdgeColor','none');
axis equal; axis tight;
title('Resolution');
colorbar;

sgtitle('Tokyo Bay Variable Resolution - Mesh Summary');
    print(gcf, fullfile(outdir, '3regions_varres_summary.png'), '-dpng', '-r150');
    fprintf('  - Saved: 3regions_varres_summary.png\n');
end

%% Overall mesh statistics
fprintf('\nMesh Statistics:\n');
fprintf('  - Total nodes: %d\n', size(m.p,1));
fprintf('  - Total elements: %d\n', size(m.t,1));
fprintf('  - Min element quality angle: %.2f degrees\n', min(min_angle_per_tri));
fprintf('  - Mean element quality angle: %.2f degrees\n', mean(min_angle_per_tri));
fprintf('  - Elements with angle < 30 degrees: %d (%.1f%%)\n', ...
        sum(min_angle_per_tri < 30), 100*sum(min_angle_per_tri < 30)/length(min_angle_per_tri));

if ENABLE_PLOTTING
    fprintf('\nAll PNG files saved to outputs/PNG directory\n');
else
    fprintf('\nPlotting was disabled - no PNG files generated\n');
end
fprintf('Mesh generation complete!\n');