clear, clc

P = xlsread('D:\Master\Empirical Methods in Finance\Homework 4\DATA_HW4.xlsx');
dates = P(:,1)+ datenum('30-Dec-1899');
P = P(:,1:end);

%Computing returns for stock/bond/rf and excess returns
R.stock = diff(P(:,1))./P(1:end-1,1);
R.bond = diff(P(:,2))./P(1:end-1,2);
Rfyear = P(2:end,3)./100;
Rf = P(2:end,3)/5200;
ExcessR.stock = R.stock-Rf; 
ExcessR.bond = R.bond-Rf;

%% 2.

% 2b. Optimal portfolio weight with lambda2 and lambda10
Variables.Mu = [mean(R.stock) mean(R.bond)]';
Variables.e = ones(2,1);
Variables.Cov = cov(R.stock,R.bond);
Variables.Rf = mean(Rf);
L = [2 10];

Weights.Static = [];

%Compute the optimal weights for different levels of risk aversion
for i = 1:2
    alpha = (1/L(i))*inv(Variables.Cov)*(Variables.Mu - Variables.e*Variables.Rf);
    alpha_rf = 1-Variables.e'*alpha;
    Weights.Static = [ Weights.Static [alpha ; alpha_rf]];
end

disp('   ')
disp('Optimal Weights')
disptable(Weights.Static,'lambda=2|lambda=10','Weight in stocks|Weight in bonds|Weight in risk free rate')

clear alpha alpha_rf i L

%% 3.

% 3.a
% Perform the Kolmogorov Test
[K.stock.H, K.stock.pval, K.stock.stat, K.stock.CV] = kstest(ExcessR.stock);
[K.bond.H, K.bond.pval, K.bond.stat, K.bond.CV] = kstest(ExcessR.bond); % The returned value of H = 1 indicates that kstest reject the null hypothesis at the default 5% significance level

disptable([K.stock.H; K.stock.stat; K.stock.pval], 'Kolmogorov-Smirnov test Stocks', 'H0|T-test|p-value')
disptable([K.bond.H; K.bond.stat; K.bond.pval], 'Kolmogorov-Smirnov test Bonds', 'H0|T-test|p-value')

% Perform the Ljung-Box test
[LB.stock.H, LB.stock.pval, LB.stock.stat, LB.stock.CV] = lbqtest(ExcessR.stock,'lags',4);
[LB.bond.H, LB.bond.pval, LB.bond.stat, LB.bond.CV] = lbqtest(ExcessR.bond,'lags',4);

disptable([LB.stock.H; LB.stock.stat; LB.stock.pval], 'Ljung-Box test Stocks', 'H0|T-test|p-value')
disptable([LB.bond.H; LB.bond.stat; LB.bond.pval], 'Ljung-Box test Bonds', 'H0|T-test|p-value')

%% 3.b estimate an AR(1) model

[AR.stock.b, AR.stock.tstat] = ols(R.stock(2:end), R.stock(1:end-1), 1);
[AR.bond.b, AR.bond.tstat] = ols(R.bond(2:end), R.bond(1:end-1), 1);

Resid.stock = R.stock(2:end)-(AR.stock.b(1)+AR.stock.b(2)*R.stock(1:end-1,1));
Resid.bond = R.bond(2:end)-(AR.bond.b(1)+AR.bond.b(2)*R.bond(1:end-1,1));

disptable([AR.stock.b(1) AR.stock.tstat(1); AR.stock.b(2) AR.stock.tstat(2)], 'Parameters Stock|t-stat', 'Mu|Phi')
disptable([AR.bond.b(1) AR.bond.tstat(1); AR.bond.b(2) AR.bond.tstat(2)], 'Parameters Bond|t-stat', 'Mu|Phi')


%% 3.c Test the ARCH effect using the LM test of Engle

[Resid2.stock.b, Resid2.stock.tstat, ~, Resid2.stock.vcv, ~, Resid2.stock.R2] = ...
    ols(Resid.stock(5:end).^2,[Resid.stock(4:end-1).^2 Resid.stock(3:end-2).^2 ...
    Resid.stock(2:end-3).^2 Resid.stock(1:end-4).^2],1);

[Resid2.bond.b, Resid2.bond.tstat, ~, Resid2.bond.vcv, ~, Resid2.bond.R2] = ...
    ols(Resid.bond(5:end).^2,[Resid.bond(4:end-1).^2 Resid.bond(3:end-2).^2 ...
    Resid.bond(2:end-3).^2 Resid.bond(1:end-4).^2],1);

[T N] = size(Resid.stock(5:end).^2);
[t n] = size(Resid.bond(5:end).^2);

LM.stock = T * Resid2.stock.R2; % If LM.stock > LM.CV => reject H0 => a1 = a2 = a3 = a4 not equal to 0
LM.bond = t * Resid2.bond.R2;
LM.CV = chi2inv(0.95,4);

disptable([LM.stock; LM.bond; LM.CV], 'LM tests 5%', 'T-test stocks|T-test bonds|Critical value')


%% 3.d Estimate the model

[GARCH.Stocks.Param, likelihood_stock, GARCH.Stocks.Sigma_hat, stderrors_stock, GARCH.Stocks.Std, scores_stock, grad_stock] = garchpq(Resid.stock, 1, 1);
[GARCH.Bonds.Param, likelihood_bond, GARCH.Bonds.Sigma_hat, stderrors_bond, GARCH.Bonds.Std, scores_bond, grad_bond] = garchpq(Resid.bond, 1, 1);

GARCH.Stocks.Tstat = GARCH.Stocks.Param./(sqrt(diag(GARCH.Stocks.Std)));
GARCH.Bonds.Tstat = GARCH.Bonds.Param./(sqrt(diag(GARCH.Bonds.Std)));

Reg = [GARCH.Stocks.Param GARCH.Stocks.Tstat GARCH.Bonds.Param GARCH.Bonds.Tstat];

disp('  ')
disp('GARCH model estimation')
disptable(Reg,'Coef.Stocks|t-stat|Coef.Bonds|t-stat','Omega|Alpha|Beta')

%ATTENTION FAUT FAIRE A+B A LA MAIN!!!!!

clear likelihood stderrors scores grad Reg


% 3e.
% Volatility Forecast
N = 52;

Sigma.Stocks = GARCH.Stocks.Param(1,1) + GARCH.Stocks.Param(2,1)*Resid.stock(end)^2 +...
    GARCH.Stocks.Param(3,1)*GARCH.Stocks.Sigma_hat(end); % 1-step ahead forecast (Sigma T+1)

Sigma.Bonds = GARCH.Bonds.Param(1,1) + GARCH.Bonds.Param(2,1)*Resid.bond(end)^2 +...
    GARCH.Bonds.Param(3,1)*GARCH.Bonds.Sigma_hat(end);   % 1-step ahead forecast (Sigma T+1)

for i = 2:N
    Sigma.Stocks(i,1) = GARCH.Stocks.Param(1,1) + (GARCH.Stocks.Param(2,1)+GARCH.Stocks.Param(3,1))...
        *Sigma.Stocks(i-1,1);
    Sigma.Bonds(i,1) = GARCH.Bonds.Param(1,1) + (GARCH.Bonds.Param(2,1)+GARCH.Bonds.Param(3,1))...
        *Sigma.Bonds(i-1,1);
end


figure('Name','Forecast Stocks')
plot(Sigma.Stocks.^(1/2));xlim([1 N]); xlabel('# weeks forecasted');
ylabel('variance forecast');
hold on ; plot((GARCH.Stocks.Param(1,1)/(1-GARCH.Stocks.Param(2,1) - GARCH.Stocks.Param(3,1)))*ones(N,1),'--r');
legend('Conditional','Unconditional')

figure('Name','Forecast Bonds')
plot(sqrt(Sigma.Bonds));xlim([1 N]); xlabel('# weeks forecasted');
ylabel('variance forecast');
hold on ; plot((GARCH.Bonds.Param(1,1)/(1-GARCH.Bonds.Param(2,1) - GARCH.Bonds.Param(3,1)))*ones(N,1),'--r');
legend('Conditional','Unconditional')