function [Hall,SensOrd]=ssi_rb(data,p,q,Sets,ref_vect,setref,r,method)

%%%Function to perform SSI covariance methods for RB cases

N=size(data{1},2);
Nsets=length(Sets);
H=cell(1,Nsets);
Href=cell(1,Nsets);
Hrov=cell(1,Nsets);
Hrov_norm=cell(1,Nsets);

%The position of the reference sensors is located
for ns=1:Nsets
    ref_vect_loc=zeros(1,length(ref_vect(1,:)));
    for i=1:length(ref_vect(1,:))
        ref_vect_loc(i)=find(ref_vect(1,i)==Sets{ns});
    end

    %The position of the roving sensors is located
    rov_vect_loc=linspace(1,length(Sets{ns}),length(Sets{ns}));
    rov_vect_loc(ref_vect_loc)=[];

    alld=[ref_vect_loc rov_vect_loc];

    rloc=size(data{ns},1);
    rref=length(ref_vect(1,:));
    %%% Method UPC. Matrix construction.
    if strcmpi(method,"UPC")
        [y_p,y_m]=hankel(data{ns},p,q,ref_vect_loc,alld);
        [~,L]=qr([y_m;y_p]',"econ");
        L=L';
        H{ns}=L(rref*q+1:end,1:rref*q);

    elseif strcmpi(method,"COVR")
        cont=1;
        R = cell(1, p+q);
        for lag=0:p+q-1
            R{cont}=data(alld,1:N-lag)*data(ref_vect_loc,lag+1:end)';
            cont=cont+1;
        end
        for row=1:(p+1)
            urow=r_j*(row-1)+1;
            for col=1:q
                ucol=rref*(col-1)+1;
                H_rr(urow:urow+r_j-1,ucol:ucol+rref-1)=R{row+col-1};
            end
        end

        %%% Method COV-DATA. Matrix construction
    elseif strcmpi(method,"COVD")
        [y_p,y_m]=Hankel(data{ns},p,q,ref_vect_loc,alld);
        H{ns}=(y_p*y_m');
    end
    %==========================================================
    
    s=svd(H{ns}, 'econ');
    rcond_rect = s(end) / s(1);
    % eps is MATLAB's internal resolution (approx. 2.22e-16).
    % If it falls below this value, the matrix is singular for the CPU.
    if rcond_rect < eps
        warning(['Warning: Hankel matrix for set %d is close to singular ...' ...
            ' or badly scaled. Results may be inaccurate. RCOND = (rcond = %e).'], ns, rcond_rect);
    end
    %==========================================================
end

rref=size(ref_vect(1,:),2);

%%% Calculation of Href(i) and Hrov(i). Reference [2] EQ(4.9) (4.19) (4.20)
for i=1:Nsets
    rrov = length(Sets{i})-rref;
    rloc = rref+rrov;

    idx_ref_loc =1:rref;
    idx_rov_loc =(rref + 1):rloc;

    steps=(0:p)'*rloc;

    idx_ref = steps+idx_ref_loc;
    idx_rov = steps+idx_rov_loc;

    idx_ref = idx_ref.'; idx_ref = idx_ref(:);
    idx_rov = idx_rov.'; idx_rov = idx_rov(:);

    Href{i} = H{i}(idx_ref, :);
    Hrov{i} = H{i}(idx_rov, :);
end

% Hall matrix assembly. Reference [2] EQ (4.24)
Href_all = zeros((p+1)*rref, Nsets*q*rref);
for i=1:Nsets
    t=q*rref*(i-1)+1:q*rref*(i);
    Href_all(:,t)=Href{i};
end

% X(i) is calculated
[~,Sref,Vref]=svd(Href_all,'econ');
Xcov_all=Sref^(1/2)*Vref';

Xcov = zeros(size(Xcov_all, 1),q*rref,Nsets);
for i=1:Nsets
    t=q*rref*(i-1)+1:q*rref*(i);
    Xcov(:,:,i)=Xcov_all(:,t);
end

% Hrov normalized is calculated.Reference [2] Algorithm 4.3
for i=1:Nsets
    if i==setref
        Hrov_norm{i}=Hrov{i};
    else
        Hrov_norm{i}=Hrov{i}*pinv(Xcov(:,:,i))*Xcov(:,:,setref);
    end
end

clearvars Hrov

%==========================================================

%%%Hall matrix assembly. Reference [2] EQ (4.24)
Hall = zeros(r*(p+1), size(Href{setref}, 2));

steps_Hall =(0:p)'*r;

%Asign reference sensors
idx_ref_Hall = steps_Hall + (1:rref);
idx_ref_Hall = idx_ref_Hall.'; idx_ref_Hall = idx_ref_Hall(:);
Hall(idx_ref_Hall, :) = Href{setref};

%Asign roving sensors
curr_offset=rref;

for Ns = 1:Nsets
    rrov_Ns = length(Sets{Ns}) - rref;
    idx_rov_loc = (curr_offset + 1):(curr_offset + rrov_Ns);

    idx_rov_Hall=steps_Hall+idx_rov_loc;
    idx_rov_Hall=idx_rov_Hall.'; idx_rov_Hall=idx_rov_Hall(:);

    Hall(idx_rov_Hall, :)=Hrov_norm{Ns};
    curr_offset=curr_offset+rrov_Ns;
end

%The order in which the Hall matrix was constructed is kept to reorganize the modal forms

SensOrd=[Sets{:}];
SensOrd=[ref_vect(1,:) SensOrd];
SensOrd=unique(SensOrd,'stable');
