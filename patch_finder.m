% Introduction:
% This script extracts regular patches of particles from a table of
% particles.
%
% Example usage:
%
% patch_finder mytable.tbl -ellipsoids [-18 18 -5.3 7 7 7;0 18 -3 7 7 7] -angular_freedom 25 -sym c4 -nearest 8 -satisfied 3 -include_extra_particles 1
%
% Author: Robert Stass
% Email: robertstass@gmail.com
% Address: OPIC, Strubi, University of Oxford.
%
% Citation: if you found this function useful for your research, please
% cite the following paper to ackownlege the author:
% tbd


function patches_table = patch_finder(table, varargin)

o = pparse.output();
disp(repmat('-',[1,60]));
disp(' ');


%%%%%%%%%%%%%%%%%%%% Default Options %%%%%%%%%%%%%%%%%%%%%%%

% find the nearest x neighbours to each particle
nearest = 8;
% Define ellipsoids that neighbours should be found in (in pixel units).
% Format: each row defines an ellipsoid. Columns 1-3 are the position of
% the centre (x y z). Columns 4-6 is the semi axis (xr yr zr).
% eg; 
ellipsoids = [-19 19 -2.3 7 7 7; 0 18 -1 7 7 7];

% Define the angular freedom allowed (this is the angle of a cone placed
% along each axis of the original particle. The angle of the neighbour must
% fall within each of these cones. 
angular_freedom = 20; %degrees

% For all the neighbours defined by 'nearest'. How many of them should
% satisfy the criteria set above.
satisfied = 3; %eg 3

% If this option is set to 1. Extra particles will be included in the
% output table. 
% These particles did not satisfy the criteria but are the neighbours of
% the succesful particles that satisfied the criteria.
include_extra_particles = 1; 
%
% symmetry used. Used to duplicate ellipsoid positions and in the angle
% check.
% only 'c' rotational symmetry supported.
sym = 'c4';
%
% Plot a region before and after to see how well it found patches. This
% makes no difference to the final table.
% set to 0 to stop this
plot_region = 0;
plot_ellipsoids = 1;
marker_size = 1.5;
%
% Define the string to append to the table filename produced. 
%This defaults to a standard format if set to []
output_append = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% input_parsing %%%%%%%%%%%%%%%%%%%%%%

    p=inputParser;
    p.addParamValue('nearest',nearest, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('ellipsoids', ellipsoids, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('sym',sym, @(x) ischar(x));
    p.addParamValue('angular_freedom',angular_freedom, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('satisfied', satisfied, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('include_extra_particles',include_extra_particles, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('plot_region',plot_region, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('plot_ellipsoids', plot_ellipsoids, @(x) isnumeric(x) || ischar(x));
    p.addParamValue('output_append', output_append, @(x) isnumeric(x) || ischar(x));

    varargin_no_minus = dynamo__parse_remove_minus(varargin);
    try
        p.parse(varargin_no_minus{:});
    catch
        disp(lasterr);
        disp('Input error. Valid parameters:');
        dynamo__multicolumn(p.Parameters,4);
        disp(' ');

        return;
    end
    %convert strings to numerical values
    q=p.Results;
    fields=fieldnames(q);
    for i=1:length(fields)
        if ~strcmp(fields{i},'sym') || ~strcmp(fields{i},'output_append') 
            q.(fields{i}) = dynamo__str2num(q.(fields{i}));
        end
    end
   
    nearest = q.nearest;
    ellipsoids = q.ellipsoids;
    sym = q.sym;
    angular_freedom = q.angular_freedom;
    satisfied = q.satisfied;
    include_extra_particles = q.include_extra_particles;
    plot_region = q.plot_region;
    plot_ellipsoids = q.plot_ellipsoids;
    output_append = q.output_append;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output file
if isequal(output_append, [])
    filestring = sprintf('_n%d_s%d_a%d_x%d',nearest, satisfied, angular_freedom, include_extra_particles);
else
    filestring = output_append;
end

% Parse the symmetry
[sym_kind,sym_parameters,ok] = dynamo_zaux_parse_sym(sym);
if isequal(ok,1)
    if isequal(sym_kind,'c')
        sym_parameter = sym_parameters;
        disp(sprintf('Using rotational symmetry of %d in angle check.', sym_parameter))
    elseif isequal(sym_kind, 'none')
        sym_parameter = sym_parameters;
        disp(sprintf('Using no rotational symmetry (%d) in angle check', sym_parameter))
    else
        disp(sprintf('Symmetry operator %s not supported', sym_kind))
        return;
    end
else
    disp(sprintf('Symmetry operator %s not understood.', sym))
    return;
end




format long g
disp(sprintf(['Finding the %d closest neighbours and checking that \n', ...
    'at least %d of them satisfy the following criteria: \n', ...
    '- the neighbour is within the %d ellipsoids defined (see plot). \n', ...
    '- its angle does not deviate by more than %d degrees from the original. \n'],nearest,satisfied,size(ellipsoids,1),angular_freedom));


ellipsoids = get_symmetric_ellipsoid_positions(ellipsoids, sym_parameter);


%%%%%%%%%%%%%%%%%%%%%% Initialise the plot %%%%%%%%%%%%%%%%%%%
%Plot the ellipsoids.
n=12; %sampling for ellipse plot
if isequal(plot_ellipsoids,1)
for i = 1:size(ellipsoids,1);
    xc = ellipsoids(i,1);
    yc = ellipsoids(i,2);
    zc = ellipsoids(i,3);
    xr = ellipsoids(i,4);
    yr = ellipsoids(i,5);
    zr = ellipsoids(i,6);
    [x_values,y_values,z_values] = ellipsoid(xc,yc,zc,xr,yr,zr,n);
    surf(x_values, y_values, z_values);%, 'edgecolor','none');
    colormap([0 1 0])
    alpha(.1)
    set(gcf, 'Position', [100,100,1000,1000]);
    view(0,90)
    set(gca,'linewidth',2)
    set(gcf, 'Renderer', 'Painters');
    hold on;
end
   

xlabel('X (pixels)');
ylabel('Y (pixels)');
zlabel('Z (pixels)');
%title('Relative positions of neighbours and ellipsoids in particle space.');
title('Relative positions of neighbours.');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55





disp('Starting calculation.');


full_table = dynamo_read(table);

alltomograms = unique(full_table(:,20));
new_table = [];
distances_list = [];
extra_particles = [];
extra_table = [];
very_large_number = 300;

% This is a bit of a mess, sorry. (I'm sure there is a better way to do this without so many for loops)
for p = 1:size(alltomograms,1)
    % for each tomogram
    tomograms_row = find(full_table(:,20) == alltomograms(p)); 
    tomo = full_table(tomograms_row,:); 
    
    allregions = unique(tomo(:,21));
    for j = 1:size(allregions,1)
        % for each region (virion)
        regions_row = find(tomo(:,21) == allregions(j)); 
        reg = tomo(regions_row,:);
        
        for n = 1:size(regions_row,1) %for each tag
            if reg(n,2)*reg(n,3) ~= 0 
                x1 = reg(n,4)+reg(n,24);
                y1 = reg(n,5)+reg(n,25);
                z1 = reg(n,6)+reg(n,26);
                data = [n reg(n,1) ones(1,6*nearest)*very_large_number];
                for i = 1:size(regions_row,1); %iterate through all the other tags for each tag.
                    if reg(i,2)*reg(i,3) ~= 0
                        x2 = reg(i,4)+reg(i,24);
                        y2 = reg(i,5)+reg(i,25);
                        z2 = reg(i,6)+reg(i,26);
                        distance = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
                        if distance ~= 0
                           for k = 1:nearest %find the nearest particles (number defined by 'nearest' variable)
                               s=(k*6)-3;
                               if distance <= data(1,s+2)
                                    %shift over existing data
                                    data = [data(1,1:s-1) 0 0 0 0 0 0 data(1,s:end)];
                                    %add new data
                                    data(s:s+2) = [i reg(i,1) distance];
                                    m = dynamo_euler2matrix(reg(n,7:9));
                                    data(s+3:s+5) = m*[x2-x1; y2-y1; z2-z1];
                                    %remove last entry
                                    data = data(1,1:end-6);
                                    break
                               end
                           end
                           
                        end
                    end
                end
                % for tag 'n' the nearest particles have been found and
                % stored in an array named data.
                % Now calculate how many of these neighbours satisfy the
                % criteria.
                score = zeros(1,nearest);
                neighbours = zeros(1,nearest);
                for k = 1:nearest %for each neighbour
                    s=(k*6)-3; %helps with indexing the array "data". 
                    particle_index = data(1,s);
                    if isequal(particle_index, very_large_number)
                        continue;
                    end
                    p = is_in_ellipsoid(ellipsoids,data(s+3),data(s+4),data(s+5)); %is it within the defined ellipsoids?
                    a = are_angles_similar(reg(n,7:9),reg(particle_index,7:9),angular_freedom, sym_parameter); %Is it's angle similar?
                    if p == 1 && a == 1 %if both conditions satisfied
                        score(k) = 1; % note this was successful.
                        neighbours(k) = data(1,s+1);		
                    end  
                end
                if sum(score) >= satisfied %if there are at least 'satisfied' many neighbours that satify the criteria.
                    distances_list = vertcat(distances_list,data); %store the data in distances_list
                    new_table = vertcat(new_table, reg(n,:)); % add the particle to the new table. 
                    % also add the particles successful neighbours to a
                    % separate list (extra particles)
                    neighbours = neighbours.*score; 
                    neighbours = neighbours(neighbours~=0);
                    extra_particles = horzcat(extra_particles,neighbours);
                end
                
            end
            %next particle.
        end
        %clean up the extra particles list to avoid duplicates. 
        extra_particles = unique(extra_particles);
        disp('Region complete.');
        %Next region.
    end
end


% If the option is turned on, add the extra particles that were neighbours
% of the selected particles.
extra_particles = setdiff(extra_particles,distances_list(:,2));
if include_extra_particles == 1 
    for m = 1:length(full_table)
        if any(extra_particles==full_table(m,1)) %if the condition is satisfied. 
            new_table = vertcat(new_table,full_table(m,:)); %write row to the new table
        end           
    end
    new_table = sortrows(new_table,1);
end
    


%Display output
before = size(full_table,1);
survived = size(distances_list,1);
extras = size(extra_particles,2);
if include_extra_particles == 1
    extras_str = 'on';
else
    extras_str = 'off';
end

disp(sprintf('\n'));
disp(sprintf(['Number of particles before: %d \n', ...
    'Number of particles after: %d \n', ...
    'An extra %d particles could be selected by accepting the successful \n', ...
    'neighbours of the surving particles bringing the total to %d particles. \n', ...
    'This option was turned %s.'],before,survived,extras,survived+extras,extras_str));

% For selected particles...
% Plot the relative positions of each neighbour relative to the particle (in particle space) 
x = [];
y = [];
z = [];
for k = 1:nearest;
    s=(k*6)-3;
    x = vertcat(x,distances_list(:,s+3));
    y = vertcat(y,distances_list(:,s+4));
    z = vertcat(z,distances_list(:,s+5));
end
if isequal(plot_ellipsoids,1)
    scatter3(x,y,z, marker_size);
    hold on
end


% Plot a before and after of one region to see the effect (set which above).
if plot_region ~= 0
    allregions = unique(full_table(:,21));
    if any(allregions==plot_region)
        tr = sprintf('reg=%d',plot_region);
        figure
        dtplot(full_table,'m','sketch','tr',tr);
        dtplot(new_table,'m','sketch','color','g','tr',tr);
        plot_central_plane('y');
    else
        disp('Region not in table');
    end
end



%Save the table
append_name = sprintf('patches%s', filestring);
new_table_name = new_filename(table,append_name); %sprintf('%s_regular_patches%s.tbl',name,filestring);
dynamo_write(new_table,new_table_name);

patches_table = new_table_name;

end


function new_name = new_filename(table,append_name)
    [path,name,ext] = fileparts(table);
    if isequal(path,'')
        path = '.';
     end
     new_name = sprintf('%s/%s_%s%s',path,name,append_name,ext);     
end






function result = is_in_ellipsoid(ellipsoids,x,y,z)
    
    for i = 1:size(ellipsoids,1);
        xc = ellipsoids(i,1);
        yc = ellipsoids(i,2);
        zc = ellipsoids(i,3);
        xr = ellipsoids(i,4);
        yr = ellipsoids(i,5);
        zr = ellipsoids(i,6);
    
        if (((x-xc)^2)/xr^2)+(((y-yc)^2)/yr^2)+(((z-zc)^2)/zr^2) < 1
            %disp('particle in ellipsoid');
            result = 1;
            return;
        else
            ;%disp('particle outside ellipsoid');
        end
            
    end
    result = 0;
    return;
end







function result = are_angles_similar(eulers1, eulers2, angular_freedom, rotational_symmetry)
    m1s = get_rotational_symmetry_equivalent_matrix(eulers1, rotational_symmetry);
    m2 = dynamo_euler2matrix(eulers2);
    for i=1:size(m1s)
        m1 = squeeze(m1s(i,:,:)); 
        mrot = m2/m1;
        px = (((mrot(2,1)^2+mrot(3,1)^2)/(2*mrot(1,1)^2))<(tan(deg2rad(angular_freedom)^2)) && mrot(1,1)>0); % is the x axis within the cone
        py = (((mrot(1,2)^2+mrot(3,2)^2)/(2*mrot(2,2)^2))<(tan(deg2rad(angular_freedom)^2)) && mrot(2,2)>0);
        pz = (((mrot(1,3)^2+mrot(2,3)^2)/(2*mrot(3,3)^2))<(tan(deg2rad(angular_freedom)^2)) && mrot(3,3)>0);
        if px == 1 && py == 1 && pz == 1 
            result = 1;
            return;
        else
            result = 0;
        end
    end
    return;
end


function ms = get_rotational_symmetry_equivalent_matrix(eulers, rotational_symmetry)
    m = dynamo_euler2matrix(eulers);
    ms = ones(rotational_symmetry,3,3);
    for s=1:rotational_symmetry
        sym_rot_eulers = [(360/rotational_symmetry)*s 0 0];
        sym_rot_m = dynamo_euler2matrix(sym_rot_eulers);
        rotm = sym_rot_m*m; 
        ms(s,:,:) = rotm;
    end    
end




function sym_ellipsoids = get_symmetric_ellipsoid_positions(ellipsoids, rotational_symmetry)
    ms = get_rotational_symmetry_equivalent_matrix([0 0 0], rotational_symmetry);
    sym_ellipsoids = [];
    for i=1:rotational_symmetry
        m = squeeze(ms(i,:,:));
        positions = ellipsoids(:,1:3);
        new_pos = transpose(m*transpose(positions));
        sym_ellipsoids = vertcat(sym_ellipsoids,(horzcat(new_pos,ellipsoids(:,4:6))));
    end
end



