function [SetsInfo,Sets,Dict]=read_input_files(filename)

%Function to read the input file containing the file names of the data and specifies how the datasets should be loaded for processing

lines=readlines(filename);
lines=lines(lines~="");

%Pre-allocation of the SetsInfo struct
SetsInfo=struct('set_id',{},'filenames',{},'parameters',{},'num_sensors_per_file', {});

current_set=NaN;
current_filename="";

%Process the file line by line
for i = 1:length(lines)
    clean_line=strip(lines(i));
    elements=strsplit(clean_line);

    fv=str2double(elements(1));

    if i==1&&~isnan(fv)
        error("The first line must contain the file name of the first set")
    end

    if isnan(fv)
        current_filename=string(elements(1));
        current_set=str2double(elements(2));

        sensor_id=str2double(elements(3));
        ang1=str2double(elements(4));
        ang2=str2double(elements(5));

        is_new_file=true;
    else
        sensor_id=str2double(elements(1));
        ang1=str2double(elements(2));
        ang2=str2double(elements(3));

        is_new_file=false;
    end

    %Check that the sensor ID and the two angles are valid numbers
    if isnan(sensor_id) || isnan(ang1) || isnan(ang2)
        error("Error in line %d: The parameters (Sensor ID, alpha and beta) must be numeric.", i);
    end

    %Check if the current set already exists
    if isempty(SetsInfo)
        idx=[];
    else
        idx=find([SetsInfo.set_id]==current_set);
    end

    %If the set does not exist, initialize its fields
    if isempty(idx)
        idx = length(SetsInfo)+1;
        SetsInfo(idx).set_id=current_set;
        SetsInfo(idx).filenames=string([]);
        SetsInfo(idx).parameters={};
        SetsInfo(idx).num_sensors_per_file=[];
    end

    %Update sensor counts and file names
    if is_new_file
        %Save the new file name
        SetsInfo(idx).filenames(end+1) = current_filename;
        %Set the sensor counter per file for this new file to 1
        SetsInfo(idx).num_sensors_per_file(end+1) = 1;
    else
        %Add 1 to the counter for the last file added to this set
        SetsInfo(idx).num_sensors_per_file(end) = SetsInfo(idx).num_sensors_per_file(end)+1;
    end

    %Save the parameters as a new row in the cell array
    SetsInfo(idx).parameters(end+1,:)={sensor_id,ang1,ang2};
end

Sets=cell(1,size(SetsInfo,2));
all_par=vertcat(SetsInfo.parameters);
filtered_par=cell2mat(all_par(:, 1:3));

[Dict,~,ic]=unique(filtered_par,'rows','stable');
current_row=1;

%The Set array is built
for k = 1:length(SetsInfo)
    if isempty(SetsInfo(k).parameters)
        error("No parameters for Set %d was specified",SetsInfo.set_id)
    end
    meas_numb_set=size(SetsInfo(k).parameters, 1);
    row_range=current_row:(current_row+meas_numb_set-1);
    ids_set=ic(row_range);

    if length(ids_set)>length(unique(ids_set))
        error("In Set %d at least one configuration is repeated",SetsInfo(k).set_id)
    end

    Sets{k}=ids_set';
    current_row =current_row+meas_numb_set;
end
