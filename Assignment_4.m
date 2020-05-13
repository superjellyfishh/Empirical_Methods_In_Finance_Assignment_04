clear, clc
%% 1. Import and Process Data

VarNames = {'SP500','JPMUS','FEDFDR'};
DataTT = readtimetable('DATA_HW4.xlsx');
DataTT.Properties.DimensionNames{1} = 'Date';
DataTT.Properties.VariableNames = VarNames;

%% 1. Compute simple returns and 'weeklized' the annualized interest rate

PricesM = DataTT{:,1:2};
ReturnsM = tick2ret(PricesM,'Method','Simple');
YieldsM = DataTT{2:end,3}./5200;
Date = DataTT.Date(2:end);

ReturnsTT = timetable(Date,ReturnsM(:,1),ReturnsM(:,2),YieldsM);
ReturnsTT.Properties.VariableNames = VarNames;

clearvars PricesM Date

%% 2. Static Allocation

lambdas = [2,10];
ExRetM = ReturnsM - YieldsM;

RowNames = {'$\mu_p$','$\sigma_p$','$\alpha_{SP500}$','$\alpha_{JPMUS}$', ...
    '$\alpha_{Rf}$'};

Results = ones(numel(RowNames),numel(lambdas));

for i = 1:numel(lambdas)
    
    lambda = lambdas(i);
    
    AlphaAs = 1/lambda*(cov(ReturnsM)\mean(ExRetM,1)');
    AlphaRf = 1 - sum(abs(AlphaAs),1);
    ReturnsP = sum(ReturnsM.*AlphaAs',2)+AlphaRf*YieldsM;
    MuP = mean(ReturnsP,1);
    SigP = std(ReturnsP);
    Results(:,i) = [MuP;SigP;AlphaAs;AlphaRf];
    
end

clearvars lambda ExRetM AlphaAs AlphaRf ReturnsP MuP SigP i

matrix2latex8(100*Results,'lambda_2_10.tex','rowlabels',RowNames, ...
    'columnlabels',{'$\lambda=2$','$\lambda=10$'},'alignment','c')

%% 3. GARCH
