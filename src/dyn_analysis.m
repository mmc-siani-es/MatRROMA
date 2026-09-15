function [clambda,W,Xi,phi,Norms,MPC,MPD,SysOrd]=dyn_analysis(Hall,max_order,sampling_rate,p,r,rref)

%%Function for calculating frequencies, damping and modal shapes for each order

%%%%Obtaining A. Reference [2] EQ (5.4) (5.5)
[Uall,Sall,Vall]=svd(Hall,'econ');
Oall=Uall*Sall^(1/2);
kappa_all=Sall^(1/2)*Vall';

%QR factorisation
Oall_up=Oall(1:p*r,1:max_order);
Oall_dw=Oall(r+1:end,1:max_order);
[Q,R]=qr(Oall_up,"econ");
S=(Q'*Oall_dw);

%The loop goes from the first pair of poles to max_order in steps
SysOrd=2:2:max_order;

clambda=cell(1,length(SysOrd));
W=cell(1,length(SysOrd));
Xi=W;
phi=W;
Norms=W;
MPC=W;
MPD=W;
dt=1/sampling_rate;
iter=0;

for i=SysOrd
    iter=iter+1;
    A=R(1:i,1:i)^-1*S(1:i,1:i);
    C=Oall(1:r,1:i);
    [psi, Val] = eig(A);
    Aq=(diag(Val))';
    Acq=log(Aq)*sampling_rate;

    %Only positive conjugated pair, excluding real poles and the negative
    %conjugated
    id_funct = find(imag(Acq)>0);
    if ~isempty(id_funct)
        clambda{iter}=Acq(id_funct);
        lambda=diag(Val);
        lambda=lambda(id_funct);
        num_poles=length(id_funct);

        %%%Modal properties
        Cpsi=C*psi(:,id_funct);
        W{iter}=abs(Acq(id_funct))/(2*pi);
        Xi{iter}=-(real(Acq(id_funct)))./abs(Acq(id_funct))*100;
        phi{iter}=Cpsi;

        %%%Transfer norm
        G_ref=kappa_all(1:i,1:rref);
        g_ref_all=psi\G_ref;
        norm_k_vector=zeros(1,num_poles);

        for k = 1:num_poles
            phi_k=Cpsi(:, k);
            g_k=g_ref_all(id_funct(k), :);
            lambda_k=lambda(k);
            omega_k=imag(Acq(id_funct(k)));
            z=exp(1i*omega_k*dt);
            norm_k_vector(k)=(norm(phi_k)*norm(g_k))/abs(z-lambda_k);
        end
        Norms{iter}=norm_k_vector;

        %%% MPC y MPD
        MPD_temp=zeros(1,num_poles);
        MPC_temp=zeros(1,num_poles);

        for k=1:num_poles
            phi_MP=Cpsi(:,k);
            %Covariances and self variances
            Sxx=sum(real(phi_MP).^2);
            Syy=sum(imag(phi_MP).^2);
            Sxy=sum(real(phi_MP).*imag(phi_MP));

            num_mpc=2*sqrt(((Sxx-Syy)/2)^2+Sxy^2);
            den_mpc=Sxx+Syy;

            %Avoid 0 divisions
            if den_mpc==0
                MPC_temp(k)= 0;
            else
                MPC_temp(k)=num_mpc/den_mpc;
            end

            %%% MPD
            % Main phase angle using atan2 to cover all 4 quadrants
            theta = 0.5*atan2(2*Sxy,Sxx-Syy);
            V=abs(phi_MP);

            phi_shift=phi_MP*exp(-1i * theta);
            deviat=abs(angle(phi_shift));

            idx_count=deviat>(pi/2);
            deviat(idx_count)=pi-deviat(idx_count);
            if sum(V) == 0
                MPD_temp(k)=0;
            else
                MPD_temp(k)=sum(V.*deviat)/sum(V);
            end
        end

        MPD{iter}=MPD_temp*(180/pi);
        MPC{iter}=MPC_temp;
    else
        clambda{iter}=[];
        W{iter}=[];
        Xi{iter}=[];
        phi{iter}=[];
        Norms{iter}=[];
        MPC{iter}=[];
        MPD{iter}=[];
    end
end

idempty=cellfun(@isempty,W);
SysOrd(idempty)=[];
clambda(idempty)=[];
W(idempty)=[];
Xi(idempty)=[];
phi(idempty)=[];
Norms(idempty)=[];
MPC(idempty)=[];
MPD(idempty)=[];
