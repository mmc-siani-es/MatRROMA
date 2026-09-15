function [data,sampling_rate,scal_fact]=process_data(SetsInfo,sampling_rate,downsampling_rate)

%%%Function that processes files containing raw data and saves them in the
%%%format required for the identification

L=inf;
data=cell(1,size(SetsInfo,2));

%The data from the files is stored in the “data” structure
for countsets=1:size(SetsInfo,2)
    clearvars P
    for countfiles=1:length(SetsInfo(countsets).filenames)
        P_temp=importdata(SetsInfo(countsets).filenames(countfiles));

        if size(P_temp,2)~=(SetsInfo(countsets).num_sensors_per_file(countfiles))
            error("The number of sensors specified in file %s does not match the number" + ...
                " of chanels in the file.",SetsInfo(countsets).filenames(countfiles))
        end

        L=min([L size(P_temp,1)]);
        if countfiles==1
            P=P_temp';
        else
            %To be able to concatenate them, we have to cut them all the the
            %size of the smaller
            P=[P(:,1:L); P_temp(1:L,:)'];
        end
    end
    data{countsets}=P;
end

if downsampling_rate>sampling_rate
    error("The downsampling frequency cannot be higher than the sampling rate. A downsampling_rate equal to the sampling rate will be taken.")
end

quad_sum=0;
n_sampl=0;

for i=1:length(data)
    data{i}=data{i}(:,1:L);
    for j=1:size(data,1)
        data{i}(j,:)=detrend(data{i}(j,:));%Data is detrended
    end

    %Data is downsampled to the desired frequency rate
    if downsampling_rate<sampling_rate
        if rem(sampling_rate,downsampling_rate)~=0
            data{i}=resample(data{i}',downsampling_rate,sampling_rate)';
        else
            factor=sampling_rate / downsampling_rate;
            dec_data=zeros(size(data{i},1),ceil(size(data{i},2)/factor));
            for ch=1:size(data{i},1)
                dec_data(ch, :)=decimate(data{i}(ch,:),factor);
            end
            data{i}=dec_data;
        end
    end
    quad_sum=quad_sum+sum(data{i}(:).^2);
    n_sampl=n_sampl+numel(data{i});
end

sigma_glob=sqrt(quad_sum/n_sampl);
scal_fact = 1/sigma_glob;
sampling_rate=downsampling_rate;

%Data is scaled to avoid ill-conditioning

for j = 1:length(data)
    data{j}=data{j}*scal_fact;
end

close all
end