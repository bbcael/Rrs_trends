%% OBSERVATIONS

clear all; close all; clc; % clear workspace
%load rrs2deg.mat; % load dataset, whether Rrs at a given resolution, Chl, SST, or anything else
%set variable of interest equal to r;
lnd = ~isnan(prod(prod(r,3),4)); % mask out land - no "(,4)" if univariate

X = [ones(20,1) (0:19)']; % design matrix
n = 20; % number of years
d = 7; % number of dimensions

for y = 1:90; % over all latitudes -- changes with resolution
    for x = 1:180; % over all longitudes -- changes with resolution
        if lnd(y,x) == 1; % over the ocean
        ts = squeeze(r(y,x,:,1:d)); % target variable for regression
        [beta,Sigma,E,CovB,logL] = mvregress(X,ts); % get multivariate multiple regression for trend
        [Ra(y,x,:,:),P1(y,x,:,:)] = corrcoef(E(1:19,:),E(2:20,:)); % calculate autocorrelation
        rhoa(y,x) = Ra(y,x,2,1); % save autocorrelation
        q1(y,x) = P1(y,x,2,1); % significance of autocorrelation
        if q1(y,x)>=0.05; % if autocorrelation is not significant -- chosen 5% level
            b = beta(2,:); c = CovB((d+1):(2*d),(d+1):(2*d)); % get coefficients & uncertainties 
            l = sqrt(sum(b.^2)); b = b./l; % get signal & normalised signal
            Z(y,x) = l./sqrt(b*c*b'); % get SNR
            q2(y,x) = q1(y,x); % autocorrelation after one Cochrane-Orcutt iteration
            q3(y,x) = q1(y,x); % same after two iterations
        elseif q1(y,x)<=0.05; % if autocorrelation is significant -- chosen 5% level -- apply Cochrane-Orcutt
            r(y,x,2:20,:) = r(y,x,2:20,:)-rhoa(y,x).*r(y,x,1:19,:); % remove autocorrelated component
            ts = squeeze(r(y,x,:,1:d)); % get target variable
            Xa = X; Xa(2:20,2) = X(2:20,2)-rhoa(y,x).*X(1:19,2); % adjust design matrix
            [beta,Sigma,E,CovB,logL] = mvregress(Xa,ts); % get multivariate multiple regression for trend
            [Rb(y,x,:,:),P2(y,x,:,:)] = corrcoef(E(1:19,:),E(2:20,:)); % calculate autocorrelation
            rhob(y,x) = Rb(y,x,2,1); % same as above
            q2(y,x) = P2(y,x,2,1); % same as above
            if q2(y,x)>=0.05; % if not significant, same as above for q1>=0.05
                b = beta(2,:); c = CovB((d+1):(2*d),(d+1):(2*d));
                l = sqrt(sum(b.^2)); b = b./l;
                Z(y,x) = l./sqrt(b*c*b');
                q3(y,x) = q2(y,x);
            elseif q2(y,x)<=0.05;  % iterate again if still significant autocorrelation
                    r(y,x,2:20,:) = r(y,x,2:20,:)-rhob(y,x).*r(y,x,1:19,:); % all same as above
                    ts = squeeze(r(y,x,:,1:d));
                    Xa = X; Xa(2:20,2) = X(2:20,2)-rhoa(y,x).*X(1:19,2);
                    Xb = Xa; Xb(2:20,2) = Xa(2:20,2)-rhob(y,x).*Xa(1:19,2);
                    [beta,Sigma,E,CovB,logL] = mvregress(Xb,ts);
                    [Rc(y,x,:,:),P3(y,x,:,:)] = corrcoef(E(1:19,:),E(2:20,:));
                    q3(y,x) = P3(y,x,2,1);
                    b = beta(2,:); c = CovB((d+1):(2*d),(d+1):(2*d));
                    l = sqrt(sum(b.^2)); b = b./l; 
                    Z(y,x) = l./sqrt(b*c*b'); % SNR
            end
        end
        end
        y
    end
end

Q = q3; % check if q3 vs q2 vs q1 is closest to a uniform distribution of p-values, or if another iteration is needed
clearvars -EXCEPT Z; save Z_2deg_Rrs.mat;
y = 89:-2:-89; y = y(1:end-5); y = repmat(y,180,1)'; % generate latitudes to calculate areal fraction with SNR>2
y(Z==0)= NaN; % remove land pixels
nansum(cosd(y(Z>2 & abs(y)<91)))./nansum(cosd(y(abs(y)<91))) % percentage of ocean area with SNR>2

%% MODEL

clear all; close all; clc; % clear workspace
% n.b. can try with adding Cochrane-Orcutt procedure as above -- one finds as expected it makes no difference to result
% load forced data as rf: 144x90x105x13 array
lnd = ~isnan(prod(prod(rf,3),4)); % mask out land
l = [412 443 488 531 547 667 678]; % observational wavelengths 
for i = 1:144; % interpolate to observed wavelengths
    for j = 1:90;
        for k = 1:size(rf,3);
            if lnd(i,j) == 1; % over the ocean
                Rif(i,j,k,:) = interp1(400:25:700,squeeze(rf(i,j,k,:)),l);        
            end
        end
    end
    i
end
clear rf i j k;
d = 7; % dimensions
X = [ones(size(Rif,3),1) (1:size(Rif,3))']; % specify design matrix
for y = 1:144; % over all longitudes
    for x = 1:88; % over all latitudes
        if lnd(y,x) == 1; % over the ocean
        ts = squeeze(Rif(y,x,:,:)); % get target variable
        if sum(isnan(ts(:)))<1060;
        [beta,Sigma,E,CovB,logL] = mvregress(X,ts); % get multivariate multiple regression for trend
        b = beta(2,:); c = CovB((d+1):(2*d),(d+1):(2*d)); % get coefficients & uncertainties
        B(y,x,:) = beta(2,:)./beta(1,:); % get slope/intercept
        sgnl(y,x) = (size(Rif,3)-1).*sqrt(sum(B(y,x,:).^2)); % get relative change over full time period
        end
        end
    end
    [y 2]
end

clearvars -EXCEPT sgnl l lnd;

% load control data as rc: 144x90x105x13 array
for i = 1:144; % interpolate to observed wavelengths
    for j = 1:90;
        for k = 1:size(rc,3);
            if lnd(i,j) == 1; % over the ocean
                Ric(i,j,k,:) = interp1(400:25:700,squeeze(rc(i,j,k,:)),l);        
            end
        end
    end
    i
end

for y = 1:144; % over all longitudes
    for x = 1:88; % over all latitudes
        if lnd(y,x) == 1; % over the ocean
        ts = squeeze(Ric(y,x,:,:)); % target variable
        [beta,Sigma,E,CovB,logL] = mvregress(X,ts); % get multivariate multiple regression for trend
        C2 = CovB(1:size(Ric,4),1:size(Ric,4)); % uncertainty in intercept
        b = squeeze(B(y,x,:))./sqrt(sum(B(y,x,:).^2)); % get normalised vector of forced trend
        U(y,x) = sqrt(b'*C2*b); % get uncertainty in intercept -- natural variability -- projected along forced trend vector
        end
    end
    [y 4]
end

ToE = 210./(sgnl'./U'); % calculate time of emergence - note 105 added to numerator because 105 added to sgnl