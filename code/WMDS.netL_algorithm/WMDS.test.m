%%
%   Input:
%         cancer: one of 14 cancer names, for exmaple 'BLCA'. Then load the
%         cancer normal and tumor matrix data, and the cancer differential
%         expression genes.
%          
%   Output:
%         a4:the drivers of this network
%

cancer = 'BLCA';
load('mRNA_lnc_priori_network');
LNC_name = importdata('../../data/WMDS.net_Run/LNC_name.txt');
filename=strcat('../../data/WMDS.net_Run',cancer,'normal.txt');
normal=importdata(filename);
filename=strcat('../../data/WMDS.net_Run',cancer,'tumor.txt');
tumor=importdata(filename);
filename=strcat('../../data/WMDS.net_Run/diff_gene_FDR',cancer,'0.01.csv');
diff_gene=importdata(filename);


text=normal.textdata(2:end,1);
for w=1:length(text)
    text{w}=text{w}(1:15);
end

IN_all=[mRNA_lnc_id_pairs(:,1);mRNA_lnc_id_pairs(:,2)];
 
IN_all=unique(IN_all);

[x1,y1]=ismember(IN_all,text);
y1(y1==0)=[];
tumor_data=tumor.data;
tumor_text=text;
IN_all_has_express=tumor_text(y1);
IN_all_tumor_expression=tumor_data(y1,:);
normal_data=normal.data;
IN_all_normal_expression=normal_data(y1,:);
IN_all_expression=[IN_all_normal_expression IN_all_tumor_expression];
[l1,l2]=find(IN_all_expression'>0);
IN_all_table=tabulate(l2);
l3=find(IN_all_table(:,2)>size(IN_all_expression,2)*0.8);
IN_all_final=IN_all_has_express(l3);

[x1,y1]=ismember(IN_all_final,text);
y1(y1==0)=[];
tumor_data=tumor.data;
tumor_text=text;
tumor_text=tumor_text(y1);
tumor_data=tumor_data(y1,:);
normal_data=normal.data;
normal_data=normal_data(y1,:);
lncRNA_text_final=tumor_text;

[~,z1]=ismember(mRNA_lnc_id_pairs(:,1),lncRNA_text_final);
[~,z2]=ismember(mRNA_lnc_id_pairs(:,2),lncRNA_text_final);
y=z1.*z2;
z=[z1 z2];
z(y==0,:)=[];
N1=length(lncRNA_text_final);
[N2,~]=size(z);
 Net=zeros(N1,N1);
for i=1:N2
         Net(z(i,2),z(i,1))=1;  %undirected gene-gene interaction network
         Net(z(i,1),z(i,2))=1;    
end
[R0,p]=xianzhu_fisher(tumor_data,normal_data);
  p(p==0)=eps;

 P6=p;
 P6(P6>=0.05)=0;
P6(P6~=0)=1; 


diff_gene_name=diff_gene.textdata(2:end,1);

N1=length(lncRNA_text_final);
[z1,~]=ismember(lncRNA_text_final,diff_gene_name);
 diff_Net=zeros(N1,N1);
 
for i=1:N1
    %i
         for j=1:N1
             diff_Net(i,j)=z1(i)+z1(j);
         end
end
diff_Net(diff_Net~=0)=1; 
         
  C=P6.*Net;   
  C=C.*diff_Net;
C=C-diag(diag(C));
 C1=C+eye(N1);
 sum_w=sum(C,2);
             sum_weight=1./sum_w;
            sum_weight(isinf(sum_weight))=10;
            w=sum_weight.^(0.01);
 f = w;
 intcon = 1:N1;
 A = -C1;
 b = -ones(N1,1);
 lb = zeros(N1,1);
 ub = ones(N1,1);
 options = optimoptions('intlinprog','Display','off');
x = intlinprog(f,intcon,A,b,[],[],lb,ub,options);

[x4,x5]=find(x); 
    [o1,o2]=find(C);
    T=tabulate(o2);
    
    [l1,l2]=sort(T(:,3),'descend') ; 
for i=1:50
    if ismember(l2(i),x4)==1
        break
    end
end
hub=l2(i);
 hub1=hub;
 for j=1:12
  hub=hub1;
 for i=1:length(hub)
 [h1,h2]=find(C(hub(i),:)==1);
 hub1=union(h2,hub1);
 end
 if length(hub1)==length(hub)
    break
end
 end
x7=intersect(x4,hub);
a4=lncRNA_text_final(x7);
filename=strcat('../../output/WMDS_output/',cancer,'_driver_lnc_singlediffexp_FDR001');
LNC_final=intersect(a4,LNC_name)
save(filename,'LNC_final');
filename=strcat('../../output/WMDS_output/',cancer,'_driver_lnc_singlediffexp_FDR001');
writetable(cell2table(LNC_final),filename)