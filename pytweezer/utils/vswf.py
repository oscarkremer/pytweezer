from .vsh import vsh

def vswf(n, m, kr, theta, phi, htype=None):
    '''
    % VSWF vector spherical wavefunctions: M_k, N_k.
    %
    % [M1,N1,M2,N2,M3,N3] = VSWF(n,m,kr,theta,phi) calculates the
    % outgoing M1,N1, incomming M2,N2 and regular M3,N3 VSWF.
    % kr, theta, phi are vectors of equal length, or scalar.
    %
    % [M,N] = VSWF(n,m,kr,theta,phi,type) calculates only the
    % requested VSWF, where type is
    %     1 -> outgoing solution - h(1)
    %     2 -> incoming solution - h(2)
    %     3 -> regular solution - j (ie RgM, RgN)
    %
    % VSWF(n, kr, theta, phi) if m is omitted, will calculate for all m.
    %
    % M,N are arrays of size length(vector_input,m) x 3
    %
    % The three components of each vector are [r,theta,phi].
    %
    % "Out of range" n and m result in return of [0 0 0]

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    #[M,N,M2,N2,M3,N3] 

    if not (isinstance(n, float) or isinstance(n, int)):
        raise TypeError('Variable \'n\' must integer of float.')
    htype = 0 if not htype else htype
    if isinstance(htype, str):        
        if htype == 'incoming':
            htype = 2
        elif htype == 'outgoing':
            htype = 1
        elif htype == 'regular':
            htype = 3
        else:
            raise ValueError('Unknown htype string')
    print(kr, theta, phi)
    B, C, P = vsh(n, m, theta, phi)

#    kr = kr(:);
#    theta = theta(:);
#    phi = phi(:);
'''

% Check the lengths
[kr,theta,phi] = matchsize(kr,theta,phi);

[B,C,P] = vsh(n,m,theta,phi);
if n > 0
    Nn = sqrt(1/(n*(n+1)));
else
    Nn = 0;
end

switch(htype)
    case 1
        if length(m)>1
            kr3=repmat(kr,[1,length(m)*3]); %makes all these suitable length
            hn=repmat(sbesselh1(n,kr),[1,length(m)*3]);
            hn1=repmat(sbesselh1(n-1,kr),[1,length(m)*3]);
        else
            kr3=threewide(kr); %makes all these suitable length
            hn=threewide(sbesselh1(n,kr));
            hn1=threewide(sbesselh1(n-1,kr));
        end
        
        M = Nn * hn .* C;
        N = Nn * ( n*(n+1)./kr3 .* hn .* P + ( hn1 - n./kr3 .* hn ) .* B );
        M2 = 0; N2 = 0; M3 = 0; N3 = 0;
    case 2
        
        if length(m)>1
            kr3=repmat(kr,[1,length(m)*3]); %makes all these suitable length
            hn=repmat(sbesselh2(n,kr),[1,length(m)*3]);
            hn1=repmat(sbesselh2(n-1,kr),[1,length(m)*3]);
        else
            kr3=threewide(kr); %makes all these suitable length
            hn=threewide(sbesselh2(n,kr));
            hn1=threewide(sbesselh2(n-1,kr));
        end
        
        M = Nn * hn .* C;
        N = Nn * ( n*(n+1)./kr3 .* hn .* P + ( hn1 - n./kr3 .* hn ) .* B );
        M2 = 0; N2 = 0; M3 = 0; N3 = 0;
    case 3
        
        if length(m)>1
            kr3=repmat(kr,[1,length(m)*3]); %makes all these suitable length
            jn=repmat(sbesselj(n,kr),[1,length(m)*3]);
            jn1=repmat(sbesselj(n-1,kr),[1,length(m)*3]);
        else
            kr3=threewide(kr); %makes all these suitable length
            jn=threewide(sbesselj(n,kr));
            jn1=threewide(sbesselj(n-1,kr));
        end
        
        M = Nn * jn .* C;
        N = Nn * ( n*(n+1)./kr3 .* jn .* P + ( jn1 - n./kr3 .* jn  ) .* B); %here is change!~!!!! get rid of jn->jn1
        M2 = 0; N2 = 0; M3 = 0; N3 = 0;
        
        if n~=1
            N(kr3==0)=0;
        else
            N(kr3==0)=2/3*Nn*( P(kr3==0) + B(kr3==0));
        end
    otherwise
        if length(m)>1
            kr3=repmat(kr,[1,length(m)*3]); %makes all these suitable length
            
            jn=repmat(sbesselj(n,kr),[1,length(m)*3]);
            jn1=repmat(sbesselj(n-1,kr),[1,length(m)*3]);
            
            hn1=repmat(sbesselh1(n,kr),[1,length(m)*3]);
            hn11=repmat(sbesselh1(n-1,kr),[1,length(m)*3]);
            
            hn2=repmat(sbesselh2(n,kr),[1,length(m)*3]);
            hn21=repmat(sbesselh2(n-1,kr),[1,length(m)*3]);
        else
            kr3=threewide(kr); %makes all these suitable length
            
            hn2=threewide(sbesselh2(n,kr));
            hn21=threewide(sbesselh2(n-1,kr));
            
            hn1=threewide(sbesselh1(n,kr));
            hn11=threewide(sbesselh1(n-1,kr));
            
            jn=threewide(sbesselj(n,kr));
            jn1=threewide(sbesselj(n-1,kr));
        end
        
        M = Nn * hn1 .* C;
        N = Nn * ( n*(n+1)./kr3 .* hn1 .* P + ( hn11 - n./kr3 .* hn1 ) .* B );
        M2 = Nn * hn2 .* C;
        N2 = Nn * ( n*(n+1)./kr3 .* hn2 .* P + ( hn21 - n./kr3 .* hn2 ) .* B );
        M3 = Nn * jn .* C;
        N3 = Nn * ( n*(n+1)./kr3 .* jn .* P + ( jn1 - n./kr3 .* jn ) .* B );
        
        if n~=1
            N3(kr3==0)=0;
        else
            N3(kr3==0)=2/3*Nn*( P(kr3==0) + B(kr3==0));
        end
end

ott.warning('external');
'''