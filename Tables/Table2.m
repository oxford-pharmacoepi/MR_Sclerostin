clear all; close all;

snps = {'rs6503468','rs9910625','rs7220711','rs66838809','rs7213935','rs80107551'};

for i = 1:length(snps)
    t_plus = readtable([snps{i} '_plus.csv']);
    t_min  = readtable([snps{i} '_min.csv']);

    % PLUS
    age = nanmean(t_plus.age);
    ageSD  = nanstd(t_plus.age);
    fem = sum(t_plus.sex == 0);
    femPC = sum(t_plus.sex == 0)/length(t_plus.sex)*100;
    mas = sum(t_plus.sex == 1);
    masPC = sum(t_plus.sex == 1)/length(t_plus.sex)*100;
    imd = nanmean(t_plus.IMD);
    imdSD = nanstd(t_plus.IMD);
    bmi = nanmean(t_plus.BMI);
    bmiSD = nanstd(t_plus.BMI);
    white = sum(ismember(t_plus.ETH,[1,1001,2001,3001,4001]));
    whitePC = white/length(t_plus.ETH)*100;
    mixed = sum(ismember(t_plus.ETH,[2,1002,2002,3002,4002]));
    mixedPC = mixed/length(t_plus.ETH)*100;
    asian = sum(ismember(t_plus.ETH,[3,1003,2003,3003,4003]));
    asianPC = asian/length(t_plus.ETH)*100;
    black = sum(ismember(t_plus.ETH,[4,1004,2004,3004,4004]));
    blackPC = black/length(t_plus.ETH)*100;
    chine = sum(t_plus.ETH == 5);
    chinePC = chine/length(t_plus.ETH)*100;
    other = sum(t_plus.ETH == 6);
    otherPC = other/length(t_plus.ETH)*100;
    plus(snps{i},age,ageSD,fem,femPC,mas,masPC,imd,imdSD,bmi,bmiSD,white,whitePC, ...
        mixed,mixedPC,asian,asianPC,black,blackPC,chine,chinePC,other,otherPC)

    % MINUS
    age = nanmean(t_min.age);
    ageSD  = nanstd(t_min.age);
    fem = sum(t_min.sex == 0);
    femPC = sum(t_min.sex == 0)/length(t_min.sex)*100;
    mas = sum(t_min.sex == 1);
    masPC = sum(t_min.sex == 1)/length(t_min.sex)*100;
    imd = nanmean(t_min.IMD);
    imdSD = nanstd(t_min.IMD);
    bmi = nanmean(t_min.BMI);
    bmiSD = nanstd(t_min.BMI);
    white = sum(ismember(t_min.ETH,[1,1001,2001,3001,4001]));
    whitePC = white/length(t_min.ETH)*100;
    mixed = sum(ismember(t_min.ETH,[2,1002,2002,3002,4002]));
    mixedPC = mixed/length(t_min.ETH)*100;
    asian = sum(ismember(t_min.ETH,[3,1003,2003,3003,4003]));
    asianPC = asian/length(t_min.ETH)*100;
    black = sum(ismember(t_min.ETH,[4,1004,2004,3004,4004]));
    blackPC = black/length(t_min.ETH)*100;
    chine = sum(t_min.ETH == 5);
    chinePC = chine/length(t_min.ETH)*100;
    other = sum(t_min.ETH == 6);
    otherPC = other/length(t_min.ETH)*100;
    min(snps{i},age,ageSD,fem,femPC,mas,masPC,imd,imdSD,bmi,bmiSD,white,whitePC, ...
        mixed,mixedPC,asian,asianPC,black,blackPC,chine,chinePC,other,otherPC)

    % SMD
    ageSMD = (nanmean(t_min.age)-nanmean(t_plus.age))/nanstd([t_min.age;t_plus.age]);
    femSMD = (nanmean(t_min.fem)-nanmean(t_plus.fem))/nanstd([t_min.fem;t_plus.fem]);
    masSMD = (nanmean(t_min.mas)-nanmean(t_plus.mas))/nanstd([t_min.mas;t_plus.mas]);
    imdSMD = (nanmean(t_min.imd)-nanmean(t_plus.imd))/nanstd([t_min.imd;t_plus.imd]);
    bmiSMD = (nanmean(t_min.bmi)-nanmean(t_plus.bmi))/nanstd([t_min.bmi;t_plus.bmi]);

end
function plus(snp,age,ageSD,fem,femPC,mas,masPC,imd,imdSD,bmi,bmiSD,white,whitePC, ...
    mixed,mixedPC,asian,asianPC,black,blackPC,chine,chinePC,other,otherPC)
writematrix(snp,'Tables.xlsx','Sheet',snp,'Range','A1')
writematrix('Reference','Tables.xlsx','Sheet',snp,'Range','A2')
writematrix('SNP (+)','Tables.xlsx','Sheet',snp,'Range','B2')
writematrix('SNP (-)','Tables.xlsx','Sheet',snp,'Range','C2')
writematrix('SMD','Tables.xlsx','Sheet',snp,'Range','D2')
writematrix('Mean (SD)','Tables.xlsx','Sheet',snp,'Range','B3')
writematrix('Mean (SD)','Tables.xlsx','Sheet',snp,'Range','C3')
writematrix('Age','Tables.xlsx','Sheet',snp,'Range','A4')
writematrix('IMD','Tables.xlsx','Sheet',snp,'Range','A5')
writematrix('BMD','Tables.xlsx','Sheet',snp,'Range','A6')
writematrix([num2str(age,"%.2f") '(' num2str(ageSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','B4')
writematrix([num2str(imd,"%.2f") '(' num2str(imdSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','B5')
writematrix([num2str(bmi,"%.2f") '(' num2str(bmiSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','B6')
writematrix('Counts (%)','Tables.xlsx','Sheet',snp,'Range','B7')
writematrix('Counts (%)','Tables.xlsx','Sheet',snp,'Range','C7')
writematrix('Sex','Tables.xlsx','Sheet',snp,'Range','A8')
writematrix('Female','Tables.xlsx','Sheet',snp,'Range','A9')
writematrix('Male','Tables.xlsx','Sheet',snp,'Range','A10')
writematrix([num2str(fem) '(' num2str(femPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B9')
writematrix([num2str(mas) '(' num2str(masPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B10')
writematrix('Ethnicity','Tables.xlsx','Sheet',snp,'Range','A11')
writematrix('White','Tables.xlsx','Sheet',snp,'Range','A12')
writematrix('Mixed','Tables.xlsx','Sheet',snp,'Range','A13')
writematrix('Asian or Asian British','Tables.xlsx','Sheet',snp,'Range','A14')
writematrix('Black or Black British','Tables.xlsx','Sheet',snp,'Range','A15')
writematrix('Chinese','Tables.xlsx','Sheet',snp,'Range','A16')
writematrix('Other ethnic groups','Tables.xlsx','Sheet',snp,'Range','A17')
writematrix([num2str(white) '(' num2str(whitePC) ')'],'Tables.xlsx','Sheet',snp,'Range','B12')
writematrix([num2str(mixed) '(' num2str(mixedPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B13')
writematrix([num2str(asian) '(' num2str(asianPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B14')
writematrix([num2str(black) '(' num2str(blackPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B15')
writematrix([num2str(chine) '(' num2str(chinePC) ')'],'Tables.xlsx','Sheet',snp,'Range','B16')
writematrix([num2str(other) '(' num2str(otherPC) ')'],'Tables.xlsx','Sheet',snp,'Range','B17')
end

function min(snp,age,ageSD,fem,femPC,mas,masPC,imd,imdSD,bmi,bmiSD,white,whitePC, ...
    mixed,mixedPC,asian,asianPC,black,blackPC,chine,chinePC,other,otherPC)
writematrix([num2str(age,"%.2f") '(' num2str(ageSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','C4')
writematrix([num2str(imd,"%.2f") '(' num2str(imdSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','C5')
writematrix([num2str(bmi,"%.2f") '(' num2str(bmiSD,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','C6')
writematrix([num2str(fem) '(' num2str(femPC,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','C9')
writematrix([num2str(mas) '(' num2str(masPC,"%.2f") ')'],'Tables.xlsx','Sheet',snp,'Range','C10')
writematrix([num2str(white) '(' num2str(whitePC) ')'],'Tables.xlsx','Sheet',snp,'Range','C12')
writematrix([num2str(mixed) '(' num2str(mixedPC) ')'],'Tables.xlsx','Sheet',snp,'Range','C13')
writematrix([num2str(asian) '(' num2str(asianPC) ')'],'Tables.xlsx','Sheet',snp,'Range','C14')
writematrix([num2str(black) '(' num2str(blackPC) ')'],'Tables.xlsx','Sheet',snp,'Range','C15')
writematrix([num2str(chine) '(' num2str(chinePC) ')'],'Tables.xlsx','Sheet',snp,'Range','C16')
writematrix([num2str(other) '(' num2str(otherPC) ')'],'Tables.xlsx','Sheet',snp,'Range','C17')
end


