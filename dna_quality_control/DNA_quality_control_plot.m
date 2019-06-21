quantity = [24.3 17.1 14.5;
12.4 17.8 7.8;
6.9 6.2 8.9;
7 7.2 11.2;
14.8 11.9 5.6;
20.5 2.7 7.7];

quality = [1.69 1.85 1.96; 
1.66 1.63 1.86;
1.78 1.87 1.96;
1.76 1.97 1.73;
1.79 1.72 2.17;
1.55 1.91 1.69];

loc = {'Neville Island'
'Point/Mon'
'Braddock'
'Point/Allegheny'
'Sharpsburg'
'Control'
};

% QUANTITY PLOT

X=1:6;
quan_mean = mean(quantity,2);
quan_std = std(quantity,[],2);

figure;

subplot(1,2,1)

bar(quan_mean)
hold on
e = errorbar(X,quan_mean,quan_std,'.');
set(gca,'XTick',1:6,'XTickLabel',loc)
grid on
title('Fig 1: DNA quantity yields for various locations and the negative control')
xlabel('Location') 
ylabel('DNA Quantity (ng/uL)') 

e.MarkerSize = 10;
% e.Color = 'red';
e.LineWidth = 1.0;



% QUALITY PLOT

qual_mean = mean(quality,2);
qual_std = std(quality,[],2);

subplot(1,2,2)

bar(qual_mean)
hold on
e = errorbar(X,qual_mean,qual_std,'.');
set(gca,'XTick',1:6,'XTickLabel',loc)
grid on
title('Fig 2: DNA quality yields for various locations and the negative control')
xlabel('Location') 
ylabel('A260/A280') 

mu = 1.8;
hline = refline([0 mu]);
hline.Color = 'r';

e.MarkerSize = 10;
% e.Color = 'red';
e.LineWidth = 1.0;


