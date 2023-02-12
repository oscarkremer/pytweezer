import numpy as np



def paraxial_transformation_matrix(paraxial_order, basis_in, basis_out, normal_mode=0):
    pass

'''

    order = int(100*normal_mode+10*basis_in[0] + basis_out[0])
    if order == 0:
        modeweights = genLG2IG(paraxial_order, 0)
    elif order == 1:
        modeweights = genLG2HG(paraxial_order)
    elif order == 2:
        modeweights = genLG2IG(paraxial_order, basis_out[1])
    elif order == 10:
        modeweights = genLG2IG(paraxial_order, 0)*genLG2vHG(paraxial_order).T
    elif order == 11:
        modeweights = genLG2HG(paraxial_order)*genLG2vHG(paraxial_order).T
    elif order == 12:
        modeweights = genLG2IG(paraxial_order, basis_out[1])*genLG2vHG(paraxial_order).T
    elif order == 20:
        modeweights = genLG2IG(paraxial_order, 0)*genLG2vIG(paraxial_order, basis_in[1]).T
    elif order == 21:
        modeweights = genLG2HG(paraxial_order)*genLG2vIG(paraxial_order,basis_in[1]).T
    elif order == 22:
        modeweights = genLG2IG(paraxial_order,basis_out[1])*genLG2vIG(paraxial_order,basis_in[1]).T
    elif order == 100:
        modeweights = np.eye(paraxial_order+1)
    elif order == 101:
        modeweights = genLG2HG(paraxial_order)*genLG2IG(paraxial_order,0).T
    elif order == 102:
        modeweights = genLG2IG(paraxial_order,basis_out[1])*genLG2IG(paraxial_order,0).T
    elif order == 110:
        modeweights = genLG2IG(paraxial_order,0)*genLG2HG(paraxial_order).T
    elif order == 111:
        modeweights = genLG2HG(paraxial_order)*genLG2HG(paraxial_order).T
    elif order == 112:
        modeweights = genLG2IG(paraxial_order,basis_out[1])*genLG2HG(paraxial_order).T
    elif order == 120:
        modeweights = genLG2IG(paraxial_order,0)*genLG2IG(paraxial_order,basis_in[1]).T
    elif order == 121:
        modeweights = genLG2HG(paraxial_order)*genLG2IG(paraxial_order,basis_in[1]).T
    elif order == 122:
        modeweights = genLG2IG(paraxial_order,basis_out[1])*genLG2IG(paraxial_order,basis_in[1]).T
    else:
        raise ValueError('Unknown order')
    if basis_out[0] == 0:
        i3_out=[]
        i2_out=[-paraxial_order:2:paraxial_order].T
        i1_out=floor((paraxial_order-abs(i2_out))/2)
    elif basis_out[0] == 1:
        i3_out=[]
        i2_out=[0:paraxial_order].T
        i1_out=paraxial_order-i2_out
    elif basis_out[0] == 2:
        i3_out=[0:paraxial_order].T
        i2_out=paraxial_order-2*(floor(i3_out/2))
        i1_out=paraxial_order*ones(paraxial_order+1,1)
        i3_out=ott.utils.iseven(i3_out)
    row_modes=[i1_out,i2_out,i3_out]
    if basis_in[0] == 0:
        i3_in=[];
        i2_in=[-paraxial_order:2:paraxial_order].T
        i1_in=floor((paraxial_order-abs(i2_in))/2)
    if basis_in[0] == 1:
        i3_in=[];
        i2_in=[0:paraxial_order].T
        i1_in=paraxial_order-i2_in
    if basis_in[0] == 2:
        i3_in=[0:paraxial_order].T
        i2_in=paraxial_order-2*(floor(i3_in/2))
        i1_in=paraxial_order*ones(paraxial_order+1,1)
        i3_in=ott.utils.iseven(i3_in)

    col_modes=[i1_in,i2_in,i3_in]
    return modeweights, col_modes, row_modes
'''
'''
def genLG2HG(order_paraxial):
    n=[0:floor(order_paraxial/2)];
    k=n;
    [N,K]=meshgrid(n,k);
    M=order_paraxial-N;
    P=min(N,M);

    normalisation_matrix=sqrt(factorial(order_paraxial-K).*factorial(K)./factorial(N)./factorial(M)/2^order_paraxial);

    summed_matrix=zeros(size(normalisation_matrix));

    summed_matrix(:,1)=factorial(M(:,1))./factorial(M(:,1)-K(:,1));
    summed_matrix(1,:)=ones([1,floor(order_paraxial/2)+1]);

    s_i=zeros(length(k)-1,floor(order_paraxial/2)+1);

    for jj=1:length(n)-1
        mm=order_paraxial-jj;
        for kk=1:length(k)-1
            s_i(kk,1)=factorial(jj)*factorial(mm)/factorial(mm-kk)/factorial(jj)*nchoosek(kk,0);            
            for ii=1:kk
                s_i(kk,ii+1)=-(jj-ii+1)/(mm-kk+ii)/(ii)*(kk-ii+1)*s_i(kk,ii);
        summed_matrix(2:end,jj+1)=sum(s_i,2);
    block_to_mirror=normalisation_matrix.*summed_matrix./factorial(K).*(-1).^(P);
    output=zeros(order_paraxial+1);
    output(1:floor(order_paraxial/2)+1,:) = [block_to_mirror(:,1:ceil(order_paraxial/2)),fliplr(block_to_mirror.*(-1).^(P.'))]; %
    output(end-ceil(order_paraxial/2)+1:end,:)=flipud(output(1:ceil(order_paraxial/2),:));
    output(ceil(order_paraxial/2)+1:end,2:2:end)=-output(ceil(order_paraxial/2)+1:end,2:2:end);
    output(2:2:end,:)=-output(2:2:end,:);
    output=output.*(1i).^repmat([0:order_paraxial].T,[1,order_paraxial+1]);
    [LGlookups,HGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
    return output, lg_lookups, hg_lookups


def genLG2vHG(order_paraxial):
    n=[0:floor(order_paraxial/2)]
    k=n
    N, K = np.meshgrid(n,k);
    M = order_paraxial-N;
    P = min(N,M)
    normalisation_matrix=sqrt(factorial(order_paraxial-K).*factorial(K)./factorial(N)./factorial(M)/2^order_paraxial);

    summed_matrix=zeros(size(normalisation_matrix));

    summed_matrix(:,1)=factorial(M(:,1))./factorial(M(:,1)-K(:,1));
    summed_matrix(1,:)=ones([1,floor(order_paraxial/2)+1]);

    s_i=zeros(length(k)-1,floor(order_paraxial/2)+1);

    for jj=1:length(n)-1
        mm=order_paraxial-jj;
        for kk=1:length(k)-1
            s_i(kk,1)=factorial(jj)*factorial(mm)/factorial(mm-kk)/factorial(jj)*nchoosek(kk,0);
            for ii=1:kk                
                s_i(kk,ii+1)=-(jj-ii+1)/(mm-kk+ii)/(ii)*(kk-ii+1)*s_i(kk,ii);
        summed_matrix(2:end,jj+1)=sum(s_i,2);
    block_to_mirror=normalisation_matrix.*summed_matrix./factorial(K).*(-1).^(P);

    output=zeros(order_paraxial+1);
    outputt=zeros(order_paraxial+1);

    outputt(1:floor(order_paraxial/2)+1,:) = [block_to_mirror(:,1:ceil(order_paraxial/2)),fliplr(block_to_mirror.*(-1).^(P.'))]; %
    outputt(end-ceil(order_paraxial/2)+1:end,:)=flipud(outputt(1:ceil(order_paraxial/2),:));

    outputt(ceil(order_paraxial/2)+1:end,2:2:end)=-outputt(ceil(order_paraxial/2)+1:end,2:2:end);
    outputt(2:2:end,:)=-outputt(2:2:end,:);
    output(1:floor((order_paraxial+1)/2),:)=1/sqrt(2)*(outputt(1:2:floor((order_paraxial+1)/2)*2,:)+outputt(2:2:floor((order_paraxial+1)/2)*2,:));

    if ~rem(order_paraxial,2)
        output(floor((order_paraxial+1)/2)+1,:)=outputt(end,:);
    output(end-floor((order_paraxial+1)/2)+1:end,:)=fliplr(flipud(output(1:floor((order_paraxial+1)/2),:)));
    output=flipud(output); 
    LGlookups, HGlookups = np.meshgrid([0:order_paraxial],[0:order_paraxial]);
    return output,LGlookups,HGlookups


def [modeweights,LGlookups,IGlookups]=genLG2IG(order_paraxial,xi)
% genLG2IG.m --- LG->IG conversion matrix for elipticity xi. A
%		parameter of xi=0 gives the non-vortex LG modes.
%		a parameter of xi=1e100 will give pretty HG modes.
%
% Usage:
%
% [modewieghts,LGlookups,IGlookups] = genLG2IG(paraxial_order,xi)

%first create the upper block... these are the fourier coefficients...
[A_n,B_n]=ott.utils.incecoefficients(order_paraxial,xi);

%prepare the index matrices for this upper block.
p=[floor(order_paraxial/2):-1:0];
l=order_paraxial-2*p;

P=repmat(p,[length(p),1]);
L=repmat(l,[length(l),1]);

Nb=(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));
Na=(1+1*(L==0)).*(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));

NA_n=flipud(Na.*A_n);
NB_n=flipud(Nb.*B_n);

%calculate the other blocks:
if rem(order_paraxial,2)==1
    bigA=[fliplr(NA_n),NA_n];
    bigB=[fliplr(NB_n),-NB_n];
    
    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+1]).';
    
else
    bigA=[fliplr(NA_n),NA_n(:,2:end)];
    bigB=[fliplr(NB_n),-NB_n(:,2:end)];
    
    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+2]).';
    modeweights(end,:)=[];
end

for ii=1:size(modeweights,1)
    modeweights(ii,:)=modeweights(ii,:)/sqrt(sum(abs(modeweights(ii,:)).^2));
end

%create imaginary matrix:
imat=repmat((-1i).^[0:order_paraxial].',[1,order_paraxial+1]);

modeweights=modeweights.*imat;

if nargout>1
    [LGlookups,IGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
end

function [modeweights,LGlookups,IGlookups]=genLG2vIG(order_paraxial,xi)
% genLG2vIG.m --- LG->IG conversion matrix for elipticity xi. A
%		parameter of xi=0 gives the vortex LG modes.
%		a parameter of xi=1e100 will give pretty vortex HG modes.
%
% Usage:
%
% [modewieghts,LGlookups,IGlookups] = genLG2vIG(paraxial_order,xi)

%first create the upper block... these are the fourier coefficients...
[A_n,B_n]=ott.utils.incecoefficients(order_paraxial,xi);

%prepare the index matrices for this upper block.
p=[floor(order_paraxial/2):-1:0];
l=order_paraxial-2*p;

P=repmat(p,[length(p),1]);
L=repmat(l,[length(l),1]);

Nb=(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));
Na=(1+1*(L==0)).*(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));

NA_n=flipud(Na.*A_n);
NB_n=flipud(Nb.*B_n);

%calculate the other blocks:
if rem(order_paraxial,2)==1
    bigA=[fliplr(NA_n),NA_n];
    bigB=[fliplr(NB_n),-NB_n];
    
    modeweightst=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+1]).';
    
else
    bigA=[fliplr(NA_n),NA_n(:,2:end)];
    bigB=[fliplr(NB_n),-NB_n(:,2:end)];
    
    modeweightst=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+2]).';
    modeweightst(end,:)=[];
end

for ii=1:size(modeweightst,1)
    
    modeweightst(ii,:)=modeweightst(ii,:)/sqrt(sum(abs(modeweightst(ii,:)).^2));
end

modeweights=zeros(size(modeweightst));

modeweights(1:floor((order_paraxial+1)/2),:)=1/sqrt(2)*(modeweightst(1:2:floor((order_paraxial+1)/2)*2,:)+modeweightst(2:2:floor((order_paraxial+1)/2)*2,:));

if ~rem(order_paraxial,2)
    modeweights(floor((order_paraxial+1)/2)+1,:)=modeweightst(end,:);
end

modeweights(end-floor((order_paraxial+1)/2)+1:end,:)=fliplr(flipud(modeweights(1:floor((order_paraxial+1)/2),:)));

if nargout>1
    [LGlookups,IGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
end
'''