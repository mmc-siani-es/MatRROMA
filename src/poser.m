function poser(filenames)

%%%Function to merge different setups with the POSER approach (with reference
% sensors)

%Pre-allocation of Dyn_par fields

omega=cell(1,length(filenames));
xi=cell(1,length(filenames));
phi=cell(1,length(filenames));
sens=cell(1,length(filenames));
ref_sens=cell(1,length(filenames));

%Load of the identified dynamic parameters of each setup
for i=1:length(filenames)
    load(filenames(i)); %#ok<*LOAD>
    omega{i}=Dyn_par.freq; %#ok<*NODEF>
    xi{i}=Dyn_par.xi;
    sens{i}=Point_DOF_order;

    if isnan(ref)
        error("POSER can be only used")
    end

    ref_sens{i}=ref;
    phi{i}=Dyn_par.phi;
end

clearvars Dyn_par Point_DOF_order ref
%========================================================
%If was not identified with a common reference the merge can not be done
if ~isequal(ref_sens{:})
    error("POSER can be only used")
end

refsens=ref_sens{i};
clearvars ref_sens

if~all(cellfun(@(x) all(ismember(refsens, x, 'rows')), sens))
    error("POSER can be only used")
end

modes_num= cellfun('size', omega, 1);

if ~all(modes_num==modes_num(1))
    error("POSER can be only used")
end

phiref=zeros(length(refsens),size(phi{1},2));

%The merge is made for each mode
for k=1:size(omega{1},2)
    Refs_mat=zeros(length(refsens),length(filenames));
    for i=1:length(filenames)
        [~,id]=ismember(refsens,sens{i},'rows');
        Refs_mat(:,i)=phi{i}(id, k); %Locate the reference positions
    end

    %Apply SVD to extract the noise-free reference mode shape
    [U, ~, ~]=svd(Refs_mat, 'econ');
    optimal_ref=U(:,1);%U(:,1) is the main direction of the reference; the rest are noise

    %The non-reference positions of each setup are normalized with
    % respect to the noise-free reference mode shape.
    for i = 1:length(filenames)
        [~, id]=ismember(refsens,sens{i},'rows');
        phi{i}(:,k)=(phi{i}(id, k)'*optimal_ref)/(phi{i}(id, k)'*phi{i}(id, k))*phi{i}(:, k);
        phiref(:,k)=phiref(:,k)+phi{i}(id, k); %The reference positions of
        %each setup are accumulated
    end
end

%The final values for the reference position are the average of the values for each
%configuration
phiref=phiref/length(filenames);

phi_global=cat(1,phi{:});
[GDL_uniq,IA, ~]=unique(vertcat(sens{:}),'rows','stable');
phi_global=phi_global(IA,:);

[~, id_uniq]=ismember(refsens,GDL_uniq,'rows');
phi_global(id_uniq,:)=phiref;


MPD_temp=zeros(size(phi_global,2),1);
MPC_temp=zeros(size(phi_global,2),1);

for i=1:size(phi_global,2)
    % Covariances and self variances
    Sxx=sum(real(phi_global(:,i)).^2);
    Syy=sum(imag(phi_global(:,i)).^2);
    Sxy=sum(real(phi_global(:,i)).*imag(phi_global(:,i)));

    num_mpc=2*sqrt(((Sxx-Syy)/2)^2+Sxy^2);
    den_mpc=Sxx+Syy;

    % Avoid 0 divisions
    if den_mpc==0
        MPC_temp(i)=0;
    else
        MPC_temp(i)=num_mpc/den_mpc;
    end
    %%% MPD
    % Main phase angle using atan2 to cover all 4 quadrants
    theta =0.5*atan2(2*Sxy,Sxx-Syy);
    V=abs(phi_global(:,i));
    phi_shift=phi_global(:,i)*exp(-1i*theta);
    deviat = abs(angle(phi_shift));

    idx_count = deviat>(pi/2);
    deviat(idx_count)=pi-deviat(idx_count);

    if sum(V)==0
        MPD_temp(i)=0;
    else
        MPD_temp(i)=sum(V.*deviat)/sum(V);
    end
end

%The parameters of the merged modes are assigned to its corresponding field in Dyn_par
Dyn_par.f=mean(cat(1,omega{:}),1);
Dyn_par.xi=mean(cat(1,xi{:}),1);
Dyn_par.phi=phi_global;
Dyn_par.MPD=MPD_temp';
Dyn_par.MPC=MPC_temp';
Point_DOF_order=GDL_uniq;
ref=refsens;

%The file containing the dynamic properties and the layout of the measured DOFs used to display the mode shapes is saved.
saved=false;
while~saved
    [file, path]=uiputfile('*.mat', 'Save selected poles as...');
    if isequal(file, 0) || isequal(path, 0)
        answer=questdlg('Are you sure you want to exit without saving?','Warning', 'Yes', 'No', 'No');
        if strcmp(answer, 'Yes')
            disp('Exiting without saving');
            saved = true;
        end
    else
        rute=fullfile(path, file);
        save(rute,'Dyn_par','Point_DOF_order',"ref");
        saved = true;
    end
end
end
