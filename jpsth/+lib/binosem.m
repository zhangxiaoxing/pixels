function [mean,sem]=binosem(pos,total)
p=pos./total;
q=1-p;
mean=p;
sem=sqrt(p*q/total);