function [y_p_f,y_m_f]=hankel(data,p,q,ref_vect,alld)

%%Function to calculate the Hankel matrix from input signals.
N=size(data,2);
rloc=size(data,1);
rref=length(ref_vect);
Nnorm=N-p-q;

%Dimensions of matrices y+ and y-: r*(p+1)xN-p-q and q*lrefxN-p-q

%Reference [1] EQ(2)
% y_p matrix construction
indices_p=(q+1:N-p)+(0:p)';
y_p_raw=data(alld,indices_p(:));
y_p=reshape(y_p_raw, rloc*(p+1),Nnorm);

% y_m matrix construction
indices_m=(0:N-p-q-1)+(q:-1:1)';
y_m_raw=data(ref_vect, indices_m(:));
y_m=reshape(y_m_raw,rref*q,Nnorm);

% Normalization
y_m_f=y_m/sqrt(Nnorm);
y_p_f=y_p/sqrt(Nnorm);
end