function [Hall,SensOrd]=ssi_rr(data,p,q,Sets,ref_vect,pairs,r,method)

%%%Function to perform SSI covariance methods for RR cases

%Pre-allocation
H_sav=cell(1,length(Sets));
rref=size(ref_vect,2);
ord=cell(1,length(Sets));
rovsens=cell(1,length(Sets));
rovsens{pairs(1,1)}=[];


%%%The first iteration is performed to obtain the Hankel matrix of the
%reference set
[H_sav{pairs(1,1)},ord{pairs(1,1)},rcond_rect]=Hankrr(data{pairs(1,1)},p,q,Sets{pairs(1,1)},ref_vect(1,:),rref,method);
if rcond_rect < eps
    warning(['Warning: Hankel matrix for set %d is close to singular ...' ...
        ' or badly scaled. Results may be inaccurate. RCOND = (rcond = %e).'], pairs(1,1), rcond_rect);
end

%%The H_sav matrix is obtained for each set
for i=1:size(pairs,1)
    [Hrr,ord{pairs(i,2)},rcond_rect]=Hankrr(data{pairs(i,2)},p,q,Sets{pairs(i,2)},ref_vect(i,:),rref,method);
    [H_sav{pairs(i,2)},rovsens{pairs(i,2)}]=assembly(q,p,Sets{pairs(i,2)},ref_vect(i,:),H_sav{pairs(i,1)},ord{pairs(i,1)},Hrr,rref);
    if rcond_rect < eps
        warning(['Warning: Hankel matrix for set %d is close to singular ...' ...
            ' or badly scaled. Results may be inaccurate. RCOND = (rcond = %e).'], pairs(i,2), rcond_rect);
    end
end

%====================================================

%%% The Hall matrix is assembled, from H_sav. Reference [2] EQ (4.24)
Hall=zeros(r*(p+1), size(H_sav{pairs(1,1)}, 2));
steps_Hall = (0:p)'*r;

%Asign reference sensors
idx_ref_Hall=steps_Hall+(1:length(ord{pairs(1,1)}));
idx_ref_Hall=idx_ref_Hall.'; idx_ref_Hall = idx_ref_Hall(:);
Hall(idx_ref_Hall,:)=H_sav{pairs(1,1)};

%Asign roving sensors
curr_offset=length(ord{pairs(1,1)}); 

for i = 1:size(pairs,1)
    rrov_Ns = length(Sets{pairs(i,2)})-rref;
    idx_rov_loc=(curr_offset+1):(curr_offset+rrov_Ns);
    idx_rov_Hall=steps_Hall+idx_rov_loc;
    idx_rov_Hall=idx_rov_Hall.'; idx_rov_Hall=idx_rov_Hall(:);
    
    idx_sav=rref+1:length(Sets{pairs(i,2)});
    idx_rov_sav=(0:p)'*length(Sets{pairs(i,2)})+idx_sav;
    idx_rov_sav=idx_rov_sav.';idx_rov_sav=idx_rov_sav(:);
    
    Hall(idx_rov_Hall, :)=H_sav{pairs(i,2)}(idx_rov_sav,:);
    curr_offset=curr_offset+rrov_Ns;
end
%_____________________________________________________

%%%The order in which the Hall matrix was constructed is kept to
%reorganize the modal shapes
SensOrd=[ord{pairs(1,1)} [rovsens{pairs(:,2)}]];
end


%%%Function to perform the hankel matrix for each set
function [H_rr,ord,rcond_rect]=Hankrr(data,p,q,Sets,ref_vect,rref,method)

r_j=size(data,1);
N=size(data,2);

%Local position of the reference sensors in the dataset
ref_vect_loc=zeros(1,rref);
for i=1:length(ref_vect)
    ref_vect_loc(i)=find(ref_vect(i)==Sets);
end

%Local position of the roving sensors in the dataset
rov_vect_loc=linspace(1,length(Sets),length(Sets));
rov_vect_loc(ref_vect_loc)=[];

alld=[ref_vect_loc rov_vect_loc];
ord=Sets(alld);

%%% Hankel matrix for UPC method
if strcmpi(method,"UPC")
    [y_p,y_m]=Hankel(data,p,q,ref_vect_loc,alld);
    [~,L]=qr([y_m;y_p]',"econ");
    L=L';
    H_rr=L(rref*q+1:end,1:rref*q);

    %%% Hankel matrix for COV-COV method
elseif strcmpi(method,"COVR")
    cont=1;
    for lag=0:p+q-1
        R(cont).cov=data(alld,1:N-lag)*data(ref_vect_loc,lag+1:end)';
        cont=cont+1;
    end
    for row=1:(p+1)
        urow=r_j*(row-1)+1;
        for col=1:q
            ucol=rref*(col-1)+1;
            H_rr(urow:urow+r_j-1,ucol:ucol+rref-1)=R(row+col-1).cov;
        end
    end

    %%% Hankel matrix for COV-data method
elseif strcmpi(method,"COVD")
    [y_p,y_m]=hankel(data,p,q,ref_vect_loc,alld);
    H_rr(:,:)=y_p*y_m';
end

s=svd(H_rr, 'econ');
rcond_rect = s(end)/s(1);

end

%%%Function to assemble the Hall matrix from the Hall of the previous iteration and the Hloc matrix of the current iteration
function [H_sav,rovsens]=assembly(q,p,Sets,ref_vect,H_sav,ord,Hrr,rref)

%Local position of the reference sensors in the master set
ref_vect_loc=zeros(1,rref);
for i=1:length(ref_vect)
    ref_vect_loc(i)=find(ref_vect(i)==ord);
end

%Extraction of Href1
rloc=length(Sets);
steps=(0:p)'*rloc;
idx_ref = steps+ref_vect_loc;
idx_ref = idx_ref.'; idx_ref = idx_ref(:);
Href1 = H_sav(idx_ref, :);


%Definition of rovsens for next set
loc_ref=zeros(1,rref);
for i=1:length(ref_vect)
    loc_ref(i)=find(ref_vect(i)==Sets);
end
rovsens=Sets;
rovsens(:,loc_ref)=[];

%Extraction of Href2 and Hrov2

rrov=length(Sets)-rref;
rloc=rref+rrov;

idx_ref_loc=1:rref;
idx_rov_loc=(rref + 1):rloc;

steps=(0:p)'*rloc;

idx_ref=steps+idx_ref_loc;
idx_rov=steps+idx_rov_loc;

idx_ref=idx_ref.';idx_ref=idx_ref(:);
idx_rov=idx_rov.';idx_rov=idx_rov(:);

Href2=Hrr(idx_ref,:);
Hrov2=Hrr(idx_rov,:);

Href_all=[Href1 Href2];

% X(i) is computed
[~,Sref,Vref]=svd(Href_all,'econ');
Xcov_all=Sref^(1/2)*Vref';

Xcov=zeros((p+1)*rref,q*rref,2);
for i=1:2
    t=q*rref*(i-1)+1:q*rref*(i);
    Xcov(:,:,i)=Xcov_all(:,t);
end

% Hrov normalized is calculated. Reference [2] Algorithm 4.3
Hrov=Hrov2*pinv(Xcov(:,:,2))*Xcov(:,:,1);

%Hsav is constructed
H_sav=zeros(rloc*(p+1),q*rref);

idx_ref_loc=1:rref;
idx_rov_loc=(rref+1):rloc;

steps=(0:p)'*rloc;

idx_ref=steps+idx_ref_loc;
idx_rov=steps+idx_rov_loc;

idx_ref=idx_ref.'; idx_ref=idx_ref(:);
idx_rov=idx_rov.'; idx_rov=idx_rov(:);


H_sav(idx_ref,:)=Href1;
H_sav(idx_rov,:)=Hrov;

end