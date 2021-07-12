% compute the coefficient of determination (r2)

function r2 = calcR2(obs,pred)

rss = nansum((obs(:)-pred(:)).^2); % residual sum of squares
tss = nansum((obs(:)-mean(obs(:))).^2); % total sum of squares
r2 = 1-rss/tss;
