function [stable_poles,modelorder_array]=stable_poles_eval(clambda,W,Xi,phi,Norm,MPC,MPD,SysOrd,tol_f,tol_xi,xi_low_lim,xi_up_lim,tol_MAC,tol_norm,cutNorm,lim_mpc,lim_mpd)

%%%Function to evaluate the stable poles to be displayed in the stabilization diagram

%Pre-allocation
stable_poles=cell(1,length(SysOrd));
modelorder_array=stable_poles;
stable_poles{1}=zeros(1,length(W{1}));
modelorder_array{1}=SysOrd(1)*ones(1,length(W{1}));


%%%The analysis compares each order i, from order 4 up to the maximum order, 
% with the previous order (i-1), searching for poles that meet the 
% stabilization criteria selected by the user. To do this, the following 
% parameters are considered: frequency deviation, damping deviation, 
% maximum and minimum damping limits, modal shape deviation, 
% modal transfer function deviation, cutoff limit in the modal transfer 
% norm function, minimum MPC value, and maximum MPD value.
for i=2:length(SysOrd)
    limNorm=sort(Norm{i},'descend');
    limNorm=limNorm(1:min(cutNorm,length(limNorm)));
    for j=1:length(W{i})
        if ismember(Norm{i}(j),limNorm) && MPC{i}(j)>=lim_mpc && MPD{i}(j)<=lim_mpd && Xi{i}(j)>xi_low_lim && Xi{i}(j)<xi_up_lim
            [~,id]=min(abs(clambda{i}(j)-clambda{i-1}));
            df=(abs((W{i}(j)-W{i-1}(id))/W{i}(j)))*100;
            dxi= (abs((Xi{i}(j)-Xi{i-1}(id))/Xi{i}(j)))*100;
            MAC=(1-(abs(phi{i}(:,j)'*phi{i-1}(:,id))^2/((phi{i}(:,j)'*phi{i}(:,j))*(phi{i-1}(:,id)'*phi{i-1}(:,id)))))*100;
            dnorm=(abs((Norm{i}(j)-Norm{i-1}(id))/Norm{i}(j)))*100;
            if df<=tol_f && dxi<=tol_xi && MAC<=tol_MAC && dnorm<=tol_norm
                stable_poles{i}(j)=1;
            else
                stable_poles{i}(j)=0;
            end
        else
            stable_poles{i}(j)=0;
        end
    end
      modelorder_array{i}=SysOrd(i)*ones(1,j);
end