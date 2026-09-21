function main(filenames,q,sampling_rate,setref,method,max_order,varargin)

%%%Master function that preprocesses the data and calls the other functions involved in the code
%%%The code requires the input signals to be arranged in rows

%%%References
%[1] Data Merging for Multi-Setup Operational Modal Analysis with
%Data-Driven SSI (DOHLER)
%[2] Subspace-based system identification and fault detection: Algorithms for large systems and application
%to structural vibration analysis (DOHLER)
%[3]REFERENCE-BASED STOCHASTIC SUBSPACE IDENTIFICATION FOR OUTPUT-ONLY
%MODAL ANALYSIS (DE ROECK)

expected_method={'UPC','COVC','COVD'};

par = inputParser;%%Parse to select
addRequired(par,'filename');
addRequired(par,'q',@(x) isnumeric(x));
addRequired(par,'sampling_rate',@(x) isnumeric(x))
addRequired(par,'setref',@(x) isnumeric(x));
addRequired(par,'method',@(x) any(validatestring(x,expected_method)));
addRequired(par,'max_order',@(x) isnumeric(x));
addParameter(par,'l',0,@(x) isnumeric(x))
addParameter(par,'ref4sinset',[],@(x) isnumeric(x));
addParameter(par,'downsampling_rate',sampling_rate,@(x) isnumeric(x))
parse(par,filenames,q,sampling_rate,setref,method,max_order,varargin{:});%%%%%

L=par.Results.l;
downsampling_rate=par.Results.downsampling_rate;
max_order=par.Results.max_order;
%(q) is the variable that defines the number of blocks in the Hankel matrix
%(L) is the variable that defines the order of the noise present in the input signals.
%(sampling_rate) is the variable that defines the sampling rate
%(setref) is the variable that defines the set of measurements taken as a reference
%(Method) defines the SSI method to use (UPC/COVD/COVC)

%The input data is adequate
[SetsInfo,Sets,Dict]=read_input_files(filenames);
[data,sampling_rate,scal_fact]=process_data(SetsInfo,sampling_rate,downsampling_rate);

%Sets are checked for correct construction
if isempty(Sets)
    error("The Sets matrix is not well constructed")
end

%Each set must have at least two sensors
for i=1:size(Sets,1)
    if length(Sets{i})<2
        error("The Sets matrix is not well constructed")
    end
end

%Reference sensors are located
if length(Sets)>1
    %The reference set must be within the number of sets.
    if setref>length(Sets)
        error('Setref is greater than the number of Sets')
    end
    [ref_vect,pairs,r]=arrangesets(Sets,setref);
else
    %Single set case
    r=max(Sets{1});
    setref=1;
    if isempty(par.Results.ref4sinset)
        warning("By not indicating reference sensors for a single set case," +...
            " by default, all sensors will act as reference, resulting in too large Hankel matrices. " +...
            "It is recommended to define reference sensors to avoid RAM overflow.") 
        ref_vect=Sets{1};
    else
        refsinset=unique(par.Results.ref4sinset,'stable');
        if all(ismember(refsinset,Sets{1}))
            ref_vect=refsinset;
        else
            error("Some of the reference sensors are not contained in Set array.")
        end
    end
end

p=q-1+L;
rref=length(nonzeros(ref_vect(1,:)'));

%The maximum order of the system is defined. Reference [2] Ap (4.6.2)
max_num_order=min([q*rref,r*(p+1)]);

if max_order>max_num_order
    max_order=max_num_order;
end

%%%Hall is obtained

% If all the reference sensors match, it is a case of Reference-Based
if all(all(ref_vect(1,:)==ref_vect))
    ref=Dict(ref_vect(1,:),:);
    [Hall,SensOrd]=ssi_rb(data,p,q,Sets,ref_vect,setref,r,method);
    % Otherwise, it is a case of Roving-References
else
    ref=NaN;
    [Hall,SensOrd]=ssi_rr(data,p,q,Sets,ref_vect,pairs,r,method);
end
%______________________________________________________________________________________

%A and C matrices are obtained from Hall. Reference [2] EQ (5.4) (5.5)
[clambda,W,Xi,phi,Norm,MPC,MPD,SysOrd]=dyn_analysis(Hall,max_order,sampling_rate,p,r,size(ref_vect,2));

%Stable poles are selected from the stability diagram
App=select_poles(data, clambda, W, Xi, phi, Norm, MPC, MPD, SysOrd, Sets, rref, sampling_rate, ref_vect, scal_fact);

Dyn_par=App.Dyn_par;
delete(App);

%The file containing the dynamic properties and the layout of the measured DOFs used to display the mode shapes is saved.
if ~isempty(Dyn_par)
    Point_DOF_order=Dict(SensOrd,:);
    saved = false;
    while ~saved
        [file, path]=uiputfile('*.mat','Save selected poles as...');
        if isequal(file, 0)||isequal(path, 0)
            answer=questdlg('Are you sure you want to exit without saving?','Warning', 'Yes', 'No', 'No'); 
            if strcmp(answer, 'Yes')
                disp('Exiting without saving');
                saved = true;
            end
        else
            rute = fullfile(path, file);
            save(rute,'Dyn_par','Point_DOF_order',"ref");
            saved = true;
        end
    end
else
    disp('No poles saved');
end
