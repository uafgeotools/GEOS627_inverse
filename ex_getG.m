function G = ex_getG(iex,m,n)
%EX_GETG return G for several examples in Aster
%
% EXAMPLES:
%    G = ex_getG(1);           % Example 1.12
%    G = ex_getG(2,100,100);   % Exercise 1-3-a
%    G = ex_getG(2,4,4);       % Exercise 1-3-e
%    G = ex_getG(2,4,20);      % Example 4.4
%    G = ex_getG(2,20,4);
%    G = ex_getG(3,210,210);   % Example 3.2
%    G = ex_getG(4,20,20);     % Example 1.6, 3.3
%    G = ex_getG(4,100,100);   % Example 1.6, 3.3
%
% all default examples at once:
%    for iex=1:4, G = ex_getG(iex); end
%
% This function is featured in lab_svd.pdf
%

bfigure = true;

% for shaw.m
%addpath('/usr/local/matlab_toolboxes/aster/cd_5.2/Lib/');

switch iex
    case 1
        stlab = 'tomography ray tracing (Ex 1.12)';
        disp('Aster Example 1.12 (and 3.1) for tomography ray tracing');
        t = sqrt(2);
        G = [1 0 0 1 0 0 1 0 0
             0 1 0 0 1 0 0 1 0
             0 0 1 0 0 1 0 0 1
             1 1 1 0 0 0 0 0 0
             0 0 0 1 1 1 0 0 0
             0 0 0 0 0 0 1 1 1
             t 0 0 0 t 0 0 0 t
             0 0 0 0 0 0 0 0 t ];
         
    case 2
        stlab = 'vertical seismic profiling (Ex 1.3)';
        disp('Aster Example 1.3, 1.9, 4.4 for the vertical seismic profiling');
        
        % m: number of slowness intervals
        % n: number of measurement intervals
        if nargin==1
            m = 100; n = m;   % Exercise 1-3-a
        end

        G = zeros(m,n);
        if n >= m
            % case n >= m (Example 4.4)
            f = n/m;
            for i = 1:m
              G(i,1:f*i) = ones(1,f*i);
            end
        else
            % case m > n (note: this is not in Aster, as far as I know)
            % construct G in n blocks, starting from the top
            f = m/n;
            for k=1:n
               inds = (k-1)*f + [1:f];
               G(inds,k) = [1:f]';
               % fill in everything to the left with m
               for kk=(k-1):-1:1
                  G(inds,kk) = f*ones(f,1);
               end
            end
        end
        
        % multiply by the depth increment
        zmin = 0; zmax = 20;    % Exercise 1-3-a
        zran = zmax-zmin;
        dz = zran / m;
        G = G*dz;
         
    case 3
        stlab = 'deconvolution (Ex 3.2)';
        disp('Aster Example 3.2 for deconvolution of instrument response (m=n)');
        % Discretizing values for M & N (210 data points)
        if nargin==1
            n = 210;
        end
        m = n;

        t = linspace(-5,100,n+1);     % time vector
        sigi = 10;
        gmax = 3.6788;

        G = zeros(m,n);
        for i = 2:m+1
          for j = 1:n
            tp = t(j)-t(i);
            if (tp > 0)
              G(i-1,j) = 0;
            else
              G(i-1,j) = -tp*exp(tp/sigi);
            end
          end
        end
        % now divide everything by the denominator
        deltat = t(2)-t(1);
        G = G/gmax * deltat;
        
    case 4
        stlab = 'Shaw slit (Ex 1.6)';
        disp('Aster Example 1.6 and Example 3.3 for the Shaw slit diffraction problem (m=n)');
        if nargin==1
            n = 20;
        end
        m = n;
        % /usr/local/matlab_toolboxes/aster/cd_5.2/Lib/shaw.m
        [G,b,x] = shaw(n);
        
        % acquire the G, m, and d originally generated using shaw.m for m = n = 20
        %load('/usr/local/matlab_toolboxes/aster/cd_5.2/Examples/chap3/ex_3_3/shaw20.mat');
        
        % acquire the G, m, and d originally generated using shaw.m for m = n = 100
        %load('/usr/local/matlab_toolboxes/aster/cd_5.2/Examples/chap3/ex_3_3/shaw100.mat'); G = G100;
end
         
[m,n] = size(G);
disp(sprintf('G is %i x %i with cond(G) = %.3e [iex = %i]',m,n,cond(G),iex));

if bfigure
    figure; imagesc(G);
    xlabel('column index k');
    ylabel('row index i');
    axis equal, axis tight;
    %set(gca,'xtick',[1:n],'ytick',[1:m]);
    colorbar;
    title(sprintf('G matrix [%i x %i] for %s',size(G),stlab));
end

%==========================================================================
