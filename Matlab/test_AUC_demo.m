clc, clear
close all
load('Segundo.mat');
I=importdata('hls_aver.txt');
W=300;
H=250;
R=reshape(I,W,H);
R=R';
figure
imshow(mat2gray(R));
N=W*H;
targets = reshape(map, 1, N); 
outputs = reshape(mat2gray(R), 1, N);
[FPR,TPR] = myPlotROC(targets, outputs);
auc1 =-trapz(FPR,TPR);
figure,plot(FPR,TPR);
xlabel('false alarm rate');
ylabel('probability of detection');
title('ROC curves of detection algorithm');
text(0.2,0.7,sprintf('AUCROC_ =%f',auc1));
[FPR1,TPR1,treh1,auc2] = myPlot3DROC(targets, outputs);
M =trapz(treh1,FPR1);
figure,plot(treh1,FPR1);
xlabel('threshold');
ylabel('False positive rate');
title('ROC curves of detection algorithm');
text(0.2,0.7,sprintf('AUCROC_ =%f',M));
set(gca,'xscale','log')
axis([10^(-4),1,0,1])