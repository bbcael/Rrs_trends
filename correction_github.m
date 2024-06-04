%% OBSERVATIONS

clear all; close all; clc; % clear workspace
cl = 0.0455; % choose confidence level -- p = 0.0455 is equivalent to SNR = 2
load AC_Rrs_2deg.mat; % load data product
r = R; % define variable on which trend detection is performed
clearvars -EXCEPT r cl; % clear workspace
oceanmask = ~isnan(prod(prod(r,3),4)); % mask out land
n = size(r,3); % number of years
d = size(r,4); % number of wavebands
X = [ones(n,1) (0:(n-1))']; % design matrix

for y = 1:size(r,1); % over all latitudes
    for x = 1:size(r,2); % over all longitudes
        if oceanmask(y,x) == 1; % over the ocean
        ts = squeeze(r(y,x,:,1:d)); % pull target variable out: Rrs(time,waveband) at each ocean grid cell 
        [beta,~,E,CovB,~] = mvregress(X,ts); % regression. outputs = estimates, residuals, covariance
        [Ra(y,x,:,:),P1(y,x,:,:)] = corrcoef(E(1:(n-1),:),E(2:n,:)); rhoa(y,x) = Ra(y,x,2,1); % compute residuals' autocorrelation
        q1(y,x) = P1(y,x,2,1); % compute p-value of this autocorrelation
        if q1(y,x)>=cl; % if autocorrelation isn't significant
            b = beta(2,:); c = CovB(2:2:(2.*d),2:2:(2.*d)); % pull out trend estimate vector and its uncertainty
            l = sqrt(sum(b.^2)); b = b./l; % calculate magnitude and direction of trend vector
            Z(y,x) = l./sqrt(b*c*b'); % calculate signal-to-noise ratio of trend
            q2(y,x) = q1(y,x); % carry p-value of autocorrelation over to first cochrane-orcutt iteration
            q3(y,x) = q1(y,x); % carry p-value of autocorrelation over to second cochrane-orcutt iteration
        elseif q1(y,x)<=cl; % if autocorrelation is significant
            r(y,x,2:n,:) = r(y,x,2:n,:)-rhoa(y,x).*r(y,x,1:(n-1),:); % adjust Rrs a la cochrane-orcutt
            ts = squeeze(r(y,x,:,1:d)); % as above
            Xa = X; Xa(2:n,2) = X(2:n,2)-rhoa(y,x).*X(1:(n-1),2); % adjust design matrix a la cochrane-orcutt
            [beta,~,E,CovB,~] = mvregress(Xa,ts); % as above
            [Rb(y,x,:,:),P2(y,x,:,:)] = corrcoef(E(1:(n-1),:),E(2:n,:)); rhob(y,x) = Rb(y,x,2,1); % as above
            q2(y,x) = P2(y,x,2,1); % as above
            if q2(y,x)>=cl; % as above
                b = beta(2,:); c = CovB(2:2:(2.*d),2:2:(2.*d)); % as above
                l = sqrt(sum(b.^2)); b = b./l; % as above
                Z(y,x) = l./sqrt(b*c*b'); % as above
                q3(y,x) = q2(y,x); % as above
            elseif q2(y,x)<=cl; % as above
                    r(y,x,2:n,:) = r(y,x,2:n,:)-rhob(y,x).*r(y,x,1:(n-1),:); % as above
                    ts = squeeze(r(y,x,:,1:d)); % as above
                    Xa = X; Xa(2:n,2) = X(2:n,2)-rhoa(y,x).*X(1:(n-1),2);% as above
                    Xb = Xa; Xb(2:n,2) = Xa(2:n,2)-rhob(y,x).*Xa(1:(n-1),2); % as above
                    [beta,~,E,CovB,~] = mvregress(Xb,ts); % as above
                    [Rc(y,x,:,:),P3(y,x,:,:)] = corrcoef(E(1:(n-1),:),E(2:n,:)); rhoc(y,x) = Rc(y,x,2,1); % as above
                    q3(y,x) = P3(y,x,2,1); % as above
                        b = beta(2,:); c = CovB(2:2:(2.*d),2:2:(2.*d)); % as above -- for final iteration
                        l = sqrt(sum(b.^2)); b = b./l; % as above
                        Z(y,x) = l./sqrt(b*c*b'); % as above
            end
        end
        end
        y./size(r,1) % print progress
    end
end
clc; % clear printed progress
ocean_lats = find(sum(oceanmask')>0); % find latitudes with ocean
empty_Z_lats = find(sum(Z')==0); % clear empty Z rows
q1(empty_Z_lats,:) = []; % as above
q2(empty_Z_lats,:) = []; % as above
q3(empty_Z_lats,:) = []; % as above
Z(empty_Z_lats,:) = []; % as above
clear ans b beta c CovB d E l n P1 P2 P3 r Ra Rb Rc rhoa rhob rhoc ts x X Xa Xb y empty_Z_lats q1 q2; % clear temporary variables
z = zeros(size(oceanmask)); z(ocean_lats,:) = Z; Z = z; clear z; % shift z-scores accordingly
Q3 = zeros(size(oceanmask)); Q3(ocean_lats,:) = q3; q3 = Q3; clear Q3; % as above
Q = q3; Q(Z==0) = NaN; % save autocorrelation p-values & NaN-ify land zeros
Z(Z==0) = NaN; % NaN-ify land zeros
clearvars -EXCEPT Z Q cl; % clear temporary variables
y = linspace(90,-90,size(Z,1)+1); y = .5.*(y(1:end-1) + y(2:end)); y = repmat(y',1,size(Z,2)); % generate latitude matrix
detected_fraction = 100.*nansum(cosd(y(Z>-norminv(cl./2))))./nansum(cosd(y(Z>0))) % fraction of ocean with significant trend
clear ans cl Q; % clear temporary variables

%% MODEL -- treated as similarly to observations as possible

clear all; close all; clc; % clear workspace
load oceanmask.mat; % load model ocean area -- excludes land and 100% coverage sea ice
load modelruns.mat; % load model Rrs interpolated to Rrs waveband centers
warning('off','all'); % suppress warnings – lots from regressions because interpolated Rrs are redundant
cl = 0.0455; % chosen confidence level
r = Rf; r0 = Rc; % define variables on which ToE analysis is performed
clear rf rc Rf Rc; % clear temporary variables

d = size(r,4); % number of wavebands
n = size(r,3); % number of years
X = [ones(size(r,3),1) (1:size(r,3))']; % specify design matrix
X20 = [ones(20,1) (1:20)']; % specify design matrix for first 20 years

for i = 1:size(r,1); % over all longitudes
    for j = 1:size(r,2); % over all latitudes
        if oceanmask(i,j) == 1; % over the ocean
        ts = squeeze(r(i,j,:,:)); % get target variable: Rrs from forced run at each location
        [beta,~,E,~,~] = mvregress(X,ts); % multivariate multiple regression for trend & trend uncertainty
        ts20 = squeeze(r(i,j,1:20,:)); % same for first 20 yaers
        [~,~,E,~,~] = mvregress(X20,ts20); % same for first 20 years
        [Ra(i,j,:,:),P1(i,j,:,:)] = corrcoef(E(1:19,:),E(2:20,:)); rhoa(i,j) = Ra(i,j,2,1);% compute significance of 20-year autocorrelation
        q1(i,j) = P1(i,j,2,1); % compute p-value of this autocorrelation
            if q1(i,j)>=cl; % if autocorrelation isn't significant
            B(i,j,:) = beta(2,:); % trend vector
            q2(i,j) = q1(i,j); % carry p-value of autocorrelation over to first cochrane-orcutt iteration
            q3(i,j) = q1(i,j); % carry p-value of autocorrelation over to second cochrane-orcutt iteration
            elseif q1(i,j)<=cl; % if autocorrelation is significant
            r(i,j,2:n,:) = r(i,j,2:n,:)-rhoa(i,j).*r(i,j,1:(n-1),:); % adjust Rrs a la cochrane-orcutt
            ts = squeeze(r(i,j,:,1:d)); % as above
            Xa = X; Xa(2:n,2) = X(2:n,2)-rhoa(i,j).*X(1:(n-1),2); % adjust design matrix a la cochrane-orcutt
            [beta,~,E,~,~] = mvregress(Xa,ts); % as above
            ts20 = squeeze(r(i,j,1:20,:)); % as above
            [~,~,E,~,~] = mvregress(X20,ts20); % as above
            [Rb(i,j,:,:),P2(i,j,:,:)] = corrcoef(E(1:19,:),E(2:20,:)); rhob(i,j) = Rb(i,j,2,1);; % as above
            q2(i,j) = P2(i,j,2,1); % as above
                if q2(i,j)>=cl; % as above
                B(i,j,:) = beta(2,:); % as above
                q3(i,j) = q2(i,j); % as above
                q4(i,j) = q2(i,j); % as above
                elseif q2(i,j)<=cl; % as above
                r(i,j,2:n,:) = r(i,j,2:n,:)-rhob(i,j).*r(i,j,1:(n-1),:); % as above
                ts = squeeze(r(i,j,:,1:d)); % as above
                Xa = X; Xa(2:n,2) = X(2:n,2)-rhoa(i,j).*X(1:(n-1),2); % as above
                Xb = Xa; Xb(2:n,2) = Xa(2:n,2)-rhob(i,j).*Xa(1:(n-1),2); % as above
                [beta,~,E,~,~] = mvregress(Xb,ts); % as above
                ts20 = squeeze(r(i,j,1:20,:)); % as above
                [~,~,E,~,~] = mvregress(X20,ts20); % as above
                [Rc(i,j,:,:),P3(i,j,:,:)] = corrcoef(E(1:19,:),E(2:20,:)); rhoc(i,j) = Rc(i,j,2,1); % as above
                q3(i,j) = P3(i,j,2,1); % as above
                B(i,j,:) = beta(2,:); % as above
                end
            end
        end
    end
    i./size(r,1) % print progress
end
clc; % clear printed progress

clearvars -EXCEPT B d X oceanmask r0 cl n; % clear temporary variables

for i = 1:size(r0,1); % as above
    for j = 1:size(r0,2); % as above
        if oceanmask(i,j) == 1; % as above
        ts = squeeze(r0(i,j,:,:)); % as above
        [beta,Sigma,E,CovB,~] = mvregress(X,ts); % as above
        c = CovB(1:2:(2.*d),1:2:(2.*d)); % covariance matrix of background from which trend emerges
        b = squeeze(B(i,j,:))./sqrt(sum(B(i,j,:).^2)); % normalized forced trend vector
        mag = sqrt(sum(B(i,j,:).^2)); % magnitude of forced trend vector
        ToE(i,j) = -norminv(cl./2).*sqrt(b'*c*b)./mag; % compute time of emergence
        end
    end
    1+i./size(r0,1) % as above
end
clc; % as above

clearvars -EXCEPT ToE oceanmask cl;
warning('on','all'); % turn warnings back on
ToE(:,87:90) = zeros(144,4); % add oceanless bands back to ToE
ToE(ToE==0) = NaN; % NaN-ify land values for ToE
ToE = ToE'; oceanmask = oceanmask'; % rotate to rows = latitude as for observations
y = linspace(-90,90,size(ToE,1)+1); y = .5.*(y(1:end-1) + y(2:end)); y = repmat(y',1,size(ToE,2)); % generate latitude matrix
median_ToE = weightedMedian(ToE(~isnan(ToE)),cosd(y(~isnan(ToE)))) % median ToE
ToE_20year_fraction = sum(cosd(y(ToE(:)<20 & oceanmask(:)==1)))./sum(cosd(y(oceanmask(:)==1))) % fraction of ocean with ToE<20 years
clearvars -EXCEPT ToE; % clear temporary variables