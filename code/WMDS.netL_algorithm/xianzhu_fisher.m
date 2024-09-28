% This code is used to implement the SSN used in the construction of the WMDS.netL algorithm



function [ dZ,p ] = xianzhu_fisher( sample,ref )
%function:construct the SSN
%   Input:
%         sample:calculated sample
%          ref:the reference samples
%   Output:
%         adjacency_matrix:the network structure
%a example
% sample=new_T(:,1);
% ref=new_N;
 [R1,P1]=corrcoef(ref');
R1(isnan(R1))=0;

[R2,P2]=corrcoef(sample');
R2(isnan(R2))=0;
Z1=1/2*log(abs((1+R1)./(1-R1)));
Z2=1/2*log(abs((1+R2)./(1-R2)));
dZ=abs(Z2-Z1)/(((1/(length(ref(1,:))-3)+1/(length(sample(1,:))-3)))^0.5);
dZ(isnan(dZ))=0;
p=1-normcdf(abs(dZ));
clear  R1 P1 R2 P2 Z1 Z2 
end