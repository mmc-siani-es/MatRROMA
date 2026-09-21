function plot_MS(MSname,coordsname,conectname)

%%% Load Data
coords=readmatrix(coordsname); % The coordinates data must be introduced in [ID, X, Y, Z] order
conect=readmatrix(conectname,'Delimiter',' '); % This file contains the relathionship of IDs

load(MSname); %#ok<LOAD>
freq=Dyn_par.f;
phi=Dyn_par.phi;

% Phi normalization
for i=1:size(phi,2)
    phi(:,i)=phi(:,i)/max(abs(phi(:,i)));
end

% ------------------------------------
num_modes=size(phi,2);
ids_nodes=coords(:,1);
coords_xyz=coords(:,2:4);
num_nodes=size(coords_xyz,1);

%%% Format check
cols_con=size(conect,2);
if cols_con==2
    type_topo='wireframe';
elseif cols_con==3
    type_topo='triangular';
else
    error('The connections are not in the correct format.');
end

conect_index=zeros(size(conect));
for j=1:numel(conect)
    rows_ind=find(ids_nodes==conect(j));
    if isempty(rows_ind)
        error('Conectivity ID not found in coordinates file.');
    end
    conect_index(j)=rows_ind;
end

%%% Model construction
displ=zeros(num_nodes,3,num_modes);

%Each DOF is associated with its corresponding coordinate
for m=1:num_modes
    for i=1:size(phi,1)
        id_meas=Point_DOF_order(i,1);
        az=Point_DOF_order(i,2);
        el=Point_DOF_order(i,3);
        idx_node=find(ids_nodes==id_meas);

        if isempty(idx_node)
            error('ID %d for modal shapes does not exist in the coordinates',id_meas);
        end

        %Each component of the modal form is broken down and added together to convert them to global coordinates
        ux=phi(i,m)*cosd(az)*cosd(el);
        uy=phi(i,m)*sind(az)*cosd(el);
        uz=phi(i,m)*sind(el);

        displ(idx_node,1,m)=displ(idx_node,1,m)+ux;
        displ(idx_node,2,m)=displ(idx_node,2,m)+uy;
        displ(idx_node,3,m)=displ(idx_node,3,m)+uz;
    end
end

%%%Figure configuration
fig = figure('Name','Mode Shapes','Color','w','Position', [100, 100, 1000, 600]);

modes_list = cell(num_modes, 1);
modes_popup{m}= cell(num_modes, 1);
for m = 1:num_modes
    modes_list{m} = sprintf('Mode %d: %.3f Hz', m, freq(m));
    modes_popup{m} = sprintf('Mode %d', m);
end

%%%% Top panels
pnl_plot = uipanel('Parent', fig, 'Title', '3D Visualization', ...
    'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', 'w', ...
    'Position', [0.02 0.22 0.96 0.75]);
ax = axes('Parent', pnl_plot, 'Position', [0.05 0.1 0.9 0.85]);
hold(ax, 'on'); grid on; axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');

if strcmp(type_topo,'wireframe')
    num_lines=size(conect_index,1);
    X_base=[coords_xyz(conect_index(:,1),1)'; coords_xyz(conect_index(:,2),1)'; nan(1,num_lines)];
    Y_base=[coords_xyz(conect_index(:,1),2)'; coords_xyz(conect_index(:,2),2)'; nan(1,num_lines)];
    Z_base=[coords_xyz(conect_index(:,1),3)'; coords_xyz(conect_index(:,2),3)'; nan(1,num_lines)];
    plot3(ax,X_base(:),Y_base(:),Z_base(:),'Color',[0.8 0.8 0.8],'LineStyle','--');
    h_geom=plot3(ax,X_base(:),Y_base(:),Z_base(:),'b-o','LineWidth',2,'MarkerFaceColor','b');
else
    trisurf(conect_index,coords_xyz(:,1),coords_xyz(:,2),coords_xyz(:,3),...
        'Parent', ax, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h_geom = trisurf(conect_index, coords_xyz(:,1), coords_xyz(:,2), coords_xyz(:,3),...
        'Parent',ax,'FaceColor','interp','EdgeColor','k');
    colormap(ax, parula);
end
margin=max(abs(displ(:)))*5;
if margin==0
    margin=1;
end
axis(ax, [min(coords_xyz(:,1))-margin, max(coords_xyz(:,1))+margin, ...
    min(coords_xyz(:,2))-margin, max(coords_xyz(:,2))+margin, ...
    min(coords_xyz(:,3))-margin, max(coords_xyz(:,3))+margin]);

%%%% Bottom panels
panel_y=0.03;
panel_h=0.15;
panel_w=0.30;
bg_color=[0.98 0.98 0.98];

% Mode selector panel
pnl_mode=uipanel('Parent',fig,'Title','Mode shape', ...
    'FontSize',10,'FontWeight','bold','BackgroundColor',bg_color, ...
    'Position', [0.02 panel_y panel_w panel_h]);

h_popup_modes=uicontrol('Parent',pnl_mode,'Style','popupmenu', ...
    'String',modes_popup,'FontSize',10, ...
    'Units', 'normalized','Position',[0.05 0.2 0.9 0.6]);

% Scale panel
pnl_scale=uipanel('Parent',fig,'Title','Deformation scale', ...
    'FontSize',10,'FontWeight','bold','BackgroundColor',bg_color, ...
    'Position', [0.35 panel_y panel_w panel_h]);

h_slider_scale=uicontrol('Parent',pnl_scale,'Style', 'slider', ...
    'Min', 0.1,'Max',8,'Value',2,...
    'Units', 'normalized','Position',[0.05 0.3 0.9 0.4]);

% Phase velocity
pnl_velo=uipanel('Parent', fig,'Title','Animation speed', ...
    'FontSize',10,'FontWeight','bold','BackgroundColor',bg_color, ...
    'Position',[0.68 panel_y panel_w panel_h]);

h_slider_velo=uicontrol('Parent',pnl_velo,'Style','slider',...
    'Min',0.01,'Max',0.5,'Value',0.15,...
    'Units','normalized','Position',[0.05 0.3 0.9 0.4]);

%%% Plotting loop
phase=0;

while isgraphics(fig)
    t_start=tic;

    actual_mode=h_popup_modes.Value;
    visual_scale=h_slider_scale.Value;
    phase_velo=h_slider_velo.Value;

    title(ax,modes_list{actual_mode},'FontSize', 14);
    actual_disp=real(displ(:,:,actual_mode)*exp(1i*phase));
    nodes_def=coords_xyz+(actual_disp*visual_scale);

    if strcmp(type_topo,'wireframe')
        num_lines=size(conect_index,1);
        X_def=[nodes_def(conect_index(:,1),1)';nodes_def(conect_index(:,2),1)';nan(1,num_lines)];
        Y_def=[nodes_def(conect_index(:,1),2)';nodes_def(conect_index(:,2),2)';nan(1,num_lines)];
        Z_def=[nodes_def(conect_index(:,1),3)';nodes_def(conect_index(:,2),3)';nan(1,num_lines)];
        set(h_geom,'XData',X_def(:),'YData',Y_def(:),'ZData',Z_def(:));
    else
        actual_magnitude=sqrt(sum(actual_disp.^2,2));
        set(h_geom, 'Vertices', nodes_def, 'FaceVertexCData', actual_magnitude);
        max_magn=max(actual_magnitude(:));
        if max_magn==0
            max_magn=1;
        end
        clim(ax,[0,max_magn]);
    end

    phase =phase+phase_velo;
    if phase>2*pi
        phase=phase-2*pi;
    end

    drawnow limitrate;
    pause(max(0,1/30-toc(t_start)));
end
