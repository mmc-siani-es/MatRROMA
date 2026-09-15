function [handle_stables, handle_uns, handle_psd]=stab_diag(data,W,stable_array,modelorder_array,origin_axis,show_psd,show_uns,rref,Sets,sampling_rate,ref_vect)

%%Function for plotting the Singular Values of the Power Spectral Density (PSD) matrix and the stabilization diagram

%%% Obtention of the Singular values of the PSD matrix
if show_psd
    NFFT=2^nextpow2(size(data{1},2)/8);
    window=hamming(NFFT/2);
    noverlap=NFFT/4; %A 50% overlap is applied
    if length(Sets)<2
        PSD=zeros(size(data{1},1),size(data{1},1),NFFT/2 + 1);
        for i = 1:size(data{1},1)
            for j = 1:size(data{1},1)
                [PSD_temp,f] = cpsd(data{1}(i, :),data{1}(j, :),window,noverlap,NFFT,sampling_rate);
                PSD(i,j,:)=PSD_temp;
                PSD(j,i,:)=conj(PSD_temp);  %The PSD is symmetrical
            end
        end
    else
        fnds=cell(1,length(Sets));
        %In the case of a reference-based case, it is not necessary to normalize the signals because there is a sensor common to all sets of measurements and
        % to find the first singular values it is sufficient to construct only one column of the PSD, that of the reference sensor with the other sensors.
        if all(all(ref_vect(1,:)==ref_vect))
            for i=1:length(Sets)
                currset=Sets{i};
                fnd=zeros(1,length(currset));
                ord_temp=unique([ref_vect(1,:) currset],'stable');
                for j=1:length(ord_temp)
                    fnd(1,j)=find(ord_temp(j)==currset);
                end
                fnds{i}=fnd;
            end

            N_freq = NFFT/2 + 1;
            ref_PSD = zeros(rref,rref,N_freq);

            %Total roving sensors rows
            total_rov_rows = 0;
            for i=1:length(Sets)
                num_sens = length(fnds{i}(1,:));
                total_rov_rows = total_rov_rows+(num_sens-rref);
            end
            rov_PSD=zeros(total_rov_rows,rref,N_freq);

            rov_idx = 1;
            % =================================================================

            for i=1:size(Sets,1)
                num_sens=length(fnds{i});
                PSD_temp=zeros(num_sens,num_sens,N_freq);

                for j=1:num_sens
                    for k=j:num_sens
                        [PSD_temp(j,k,:),f]=cpsd(data{i}(fnds{i}(1,j),:), ...
                            data{i}(fnds{i}(1,k),:), window, noverlap, NFFT, sampling_rate);
                    end
                end

                %Constructuion of complete ref:psd matrix
                ref_PSD=ref_PSD+PSD_temp(1:rref,1:rref,:);

                %Iterative construction of rov PDS matrix for each set with
                %rov_idx
                num_rov_current = num_sens - rref;
                if num_rov_current > 0
                    rov_PSD(rov_idx:rov_idx+num_rov_current-1,1:rref,:)=PSD_temp(rref+1:end,1:rref,:);
                    rov_idx=rov_idx+num_rov_current;
                end
            end

            for i=1:N_freq
                ref_PSD(:,:,i)=ref_PSD(:,:,i)+triu(ref_PSD(:,:,i),1)';
            end

            PSD=[ref_PSD/size(Sets,1);rov_PSD];

        else
            %Roving References Case: Averaged Spectrum (ANPSD)
            % Since there is no common reference, the PSD and the
            % 1st singular value of each set are calculated independently and averaged.
            S_avg=zeros(1, NFFT/2+1);

            for s = 1:length(Sets)
                data_mat=data{s};
                n_sens=size(data_mat, 1);

                % PSD for each set
                PSD_set=zeros(n_sens,n_sens,NFFT/2+1);
                for i=1:n_sens
                    for j = i:n_sens
                        [PSD_temp, f]=cpsd(data_mat(i,:),data_mat(j,:),window,noverlap,NFFT,sampling_rate);
                        PSD_set(i,j,:)=PSD_temp;
                        if i ~= j
                            PSD_set(j,i,:)=conj(PSD_temp);
                        end
                    end
                end

                %Temporary vector to store the first singular value of this set
                S_set_curve=zeros(1, NFFT/2 + 1);

                for freq_idx = 1:NFFT/2+1
                    [~,Sb,~]=svd(PSD_set(:,:,freq_idx),'econ');
                    S_set_curve(freq_idx)=Sb(1,1);
                end

                S_set_norm=S_set_curve/max(S_set_curve);
                S_avg =S_avg+S_set_norm;

            end
            S=S_avg/size(Sets,1);
        end
    end

    if ~exist('S',"var")
        %The singular values of the PSD are obtained.
        for i=1:NFFT/2+1
            [~,Sb,~]=svd(PSD(:,:,i),'econ');
            S(:,i)=diag(Sb);
        end
    end
    S=S(1,:)./max(S(1,:)); % normalized power spectral density

    %Convert to decibels to preserve the logarithmic scale on a linear axis
    S_dB=10*log10(S(1,:));
    min_dB=min(S_dB(isfinite(S_dB)));% Ignore possible -Inf values due to zeros
    max_dB=max(S_dB); %Since S is normalized to 1, max_dB will be 0
    S_norm=(S_dB-min_dB)/(max_dB-min_dB);

    %Multiply so that the highest peak reaches 80% of the maximum order
    S_plot=S_norm*modelorder_array{length(modelorder_array)}(1,1)*0.8;
end


%%%Stable and non-stable poles are stored in two separate vectors for later plotting.

evalmat=[[W{:}]',[modelorder_array{:}]'];
stid=[stable_array{:}]'~=0;

stablepoles = evalmat(stid,:);
nonstablepoles = evalmat(~stid,:);

%==========================================
%%%Dynamic plotting

cla(origin_axis);
hold(origin_axis, 'on');

handle_psd=[];
handle_uns=[];
%Plot PSD only if requested
if show_psd
    handle_psd=semilogy(origin_axis, f, S_plot,'Color','#B0B0B0','LineWidth',0.3,'HitTest','off');
end

yyaxis(origin_axis, 'left');

%Plot unstable poles only if requested
if show_uns
    handle_uns=plot(origin_axis, nonstablepoles(:,1), nonstablepoles(:,2), 'o',...
        'Color','#666666','MarkerSize',6,'LineStyle','none','HitTest','on');
end

handle_stables=plot(origin_axis,stablepoles(:,1),stablepoles(:,2),'o', ...
    'LineWidth',2,'Color', '#228B22', 'MarkerSize',6,'LineStyle','none','HitTest','on');

hold(origin_axis,'off');

if ~isempty(handle_psd) && isgraphics(handle_psd)
    %Get all the objects drawn on the axes
    allChildren=origin_axis.Children;

    %Check out the PSD's position
    idx=(allChildren==handle_psd);

    if any(idx)
        %Extract the remaining objects (stable and unstable poles)
        otherChildren=allChildren(~idx);

        %Reassign the PSD to the background
        origin_axis.Children=[otherChildren;handle_psd];
    end
end

% Visual configuration of the axes
xlabel(origin_axis, "Frequency (Hz)",'FontSize', 13,'FontWeight','bold',"Interpreter","latex");
ylabel(origin_axis, "System Order",'FontSize',13,'FontWeight','bold',"Interpreter","latex");
xlim(origin_axis, [0,sampling_rate/2]);
ylim(origin_axis,[0,max(cellfun(@max, modelorder_array))+2]);
origin_axis.YAxis(1).Color='k';
origin_axis.YAxis(2).Color='k';
origin_axis.YAxis(2).TickLabels="";
origin_axis.YAxis(2).TickLength=[0 0];

end