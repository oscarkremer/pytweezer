import numpy as np
from .legendre_row import legendre_row
from .match_size import match_size

def spharm(n, m, theta, phi=None):
    '''
    % SPHARM scalar spherical harmonics and angular partial derivatives.
    %
    % Y = SPHARM(n,m,theta,phi) calculates scalar spherical harmonics.
    %
    % [Y,Ytheta,Yphi] = SPHARM(n,m,theta,phi) additionally, calculates
    % the angular partial derivatives dY/dtheta and 1/sin(theta)*dY/dphi.
    %
    % SPHARM(n,theta,phi) as above but for all m.
    %
    % Scalar n for the moment.
    % 
    % If scalar m is used Y is a vector of length(theta,phi) and is
    % completely compatible with previous versions of the toolbox. If vector m
    % is present the output will be a matrix with rows of length(theta,phi) for
    % m columns.
    %
    % "Out of range" n and m result in return of Y = 0

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if not (isinstance(n, float) or (isinstance(n, int)):
        raise TypeError('Input parameter \'n\' must be scalar.')
    if not phi:    
        phi = theta
        theta = m
    mi = m
    m = np.arange(-n:n+1)
    index_bigger = np.arghwere(abs(m) <=n)
    if index_bigger.size:
       m = m([tuple(index_bigger.T)])

    [theta, phi] = match_size(theta, phi)
    #input_length = 



    #pnm = legendrerow(n,theta);


    #pnm = pnm(abs(m)+1,:); %pick the m's we potentially have.

    #[phiM,mv]=meshgrid(phi,m);

    #pnm = [(-1).^mv(m<0,:).*pnm(m<0,:);pnm(m>=0,:)];

    #expphi = exp(1i*mv.*phiM);
    #Y = pnm .* expphi;
    #if nargout <= 1
    #    Y=Y.';
    #    ott.warning('external');
    #    return

    #expplus = exp(1i*phiM);
    #expminus = exp(-1i*phiM);

    '''
    ymplus=[Y(2:end,:);zeros(1,length(theta))];
    ymminus=[zeros(1,length(theta));Y(1:end-1,:)];

    Ytheta = sqrt((n-mv+1).*(n+mv))/2 .* expplus .* ymminus ...
            - sqrt((n-mv).*(n+mv+1))/2 .* expminus .* ymplus;

    Y2 = spharm(n+1,theta,phi).';

    ymplus=Y2(3:end,:);
    ymminus=Y2(1:end-2,:);

    Yphi = 1i/2 * sqrt((2*n+1)/(2*n+3)) * ...
    ( sqrt((n+mv+1).*(n+mv+2)) .* expminus .* ymplus ...
    + sqrt((n-mv+1).*(n-mv+2)) .* expplus .* ymminus );

    Y=Y(n+mi+1,:).';
    Yphi=Yphi(n+mi+1,:).';
    Ytheta=Ytheta(n+mi+1,:).';

    return Y, Ytheta, Yphi
    '''