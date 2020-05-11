clear, clc
%% 1. Import and Process Data
%% Import the data, extracting spreadsheet dates in Excel serial date format
[~, ~, raw, dates] = xlsread('DATA_HW4.xlsx','Feuil1','A3:D941','',@convertSpreadsheetExcelDates);
raw = raw(:,[2,3,4]);
dates = dates(:,1);
data = reshape([raw{:}],size(raw));
DATAHW4 = table;
% Allocate imported array to column variable names
DATAHW4.Name = datetime([dates{:,1}].', 'ConvertFrom', 'Excel');
DATAHW4.SP500 = data(:,1);
DATAHW4.FI = data(:,2);
DATAHW4.IR = data(:,3);
% Clear temporary variables
clearvars data raw dates;

%% 1. Compute simple returns and 'weeklized' the annualized interest rate

R_SP500 = tick2ret(DATAHW4.SP500);
R_bonds = tick2ret(DATAHW4.FI);
R_IR = ((DATAHW4.IR/100/52));
R_IR = R_IR(2:end); % So that we have the same dimensions

% The excess returns:
R_SP500_exc = R_SP500 - R_IR;
R_FI_exc = R_bonds - R_IR;

%% 2a.
mu_SP500 = mean(R_SP500_exc);
mu_FI = mean(R_FI_exc);

mu = [mu_SP500 mu_FI]';
sigma = cov(R_SP500, R_bonds);
e = ones(2,1);
alphas = zeros(3,2);

% Lambda = 2
lambda = 2
n = size(mu_SP500, 1);
weights2 = zeros(n + 1, 1);
weights2(1:n) = (1/lambda) * (sigma \ mu);
weights2(n + 1) = 1 - sum(weights2(1:n));

% Lambda = 10
lambda = 10
n = size(mu_SP500, 1);
weights10 = zeros(n + 1, 1);
weights10(1:n) = (1/lambda) * (sigma \ mu);
weights10(n + 1) = 1 - sum(weights10(1:n));





