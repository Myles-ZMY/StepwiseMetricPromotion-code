clc
clear
rand('state',1) 
addpath('./code/');
% addpath('./feat/');
  load(['PRID2011_LOMO_10_iteration.mat']);
% load(['PRID_GOGdim600_10_iteration.mat']);
sigma = 10000;
numK=10;
Knear = 10;

for iter=1:10
   numPCA = 600;
fea_tr_a =   (  feat{iter}.fea_tr_a  );  
fea_tr_b =    ( feat{iter}.fea_tr_b  );  
fea_te_a =    (  feat{iter}.fea_te_a ) ; 
fea_te_b =    (  feat{iter}.fea_te_b ) ;
tr_label_a = feat{iter}.tr_label_a;
tr_label_b = feat{iter}.tr_label_b;
te_label_a = feat{iter}.te_label_a;
te_label_b = feat{iter}.te_label_b;
tr_id = feat{iter}.tr_id ;
te_id = feat{iter}.te_id ;

  t0 = clock;
 [W1t, M1t,inCov1,exCov1] = XQDA(fea_tr_a', fea_tr_a',tr_label_a , tr_label_a );
 [W2t, M2t,inCov2,exCov2] = XQDA(fea_tr_b', fea_tr_b',tr_label_b', tr_label_b');
  M1 = W1t*M1t*W1t';
  
  M2 = W2t*M2t*W2t';
  


kk=1;
while (1)
dis11 = MahDist(M1, fea_tr_a', fea_tr_a')+MahDist(M2, fea_tr_a', fea_tr_a');
dis12 =  MahDist(M1, fea_tr_b', fea_tr_a')+MahDist(M2, fea_tr_b', fea_tr_a') ;
dis21 =  MahDist(M2, fea_tr_a', fea_tr_b')+MahDist(M1, fea_tr_a', fea_tr_b') ;
dis22 = MahDist(M2, fea_tr_b', fea_tr_b')+ MahDist(M1, fea_tr_b', fea_tr_b');
for ii=1:89
    for jj=1:89
%       dist_maxpooling2(jj,ii) = MahDist(M2, max(fea_te_b(:,find(te_label_b==te_id(jj))),[],2)'*W2,  max(fea_te_a(:,find(te_label_a==te_id(ii))),[],2)'*W2);
        temp11 = dis11(find(tr_label_a==tr_id(jj)),find(tr_label_a==tr_id(ii)));
        temp12 = dis12(find(tr_label_b==tr_id(jj)),find(tr_label_a==tr_id(ii)));
        temp21 = dis21(find(tr_label_a==tr_id(ii)),find(tr_label_b==tr_id(jj)));
        temp22 = dis22(find(tr_label_b==tr_id(ii)),find(tr_label_b==tr_id(jj)));
        Dis11(jj,ii)= sum(temp11(:))/length(find(tr_label_a==tr_id(jj)))/length(find(tr_label_a==tr_id(ii)));
        Dis12(jj,ii)= sum(temp12(:))/length(find(tr_label_b==tr_id(jj)))/length(find(tr_label_a==tr_id(ii)));
        Dis21(jj,ii)= sum(temp21(:))/length(find(tr_label_b==tr_id(jj)))/length(find(tr_label_a==tr_id(ii)));
        Dis22(jj,ii)= sum(temp22(:))/length(find(tr_label_b==tr_id(jj)))/length(find(tr_label_b==tr_id(ii)));
    end
end

ptmp12 = zeros(89,89);ptmp21 = zeros(89,89);
ptmp11 = zeros(89,89);ptmp22 = zeros(89,89);
for icol = 1:89 
 [val12 pos12 ] =  sort(Dis12(:,icol));
  ptmp12(pos12(1:Knear),icol) = 1;
 [ val21 pos21] = sort(Dis21(:,icol));
  ptmp21(pos21(1:Knear),icol) = 1;
   [val11 pos11 ] =  sort(Dis11(:,icol));
  ptmp11(pos11(1:numK),icol) = 1;
 [ val22 pos22] = sort(Dis22(:,icol));
  ptmp22(pos22(1:numK),icol) = 1; 
end



ptmp = ptmp12.*ptmp21';
% % %
% ppp = ptmp';
s12 = exp(-abs(Dis12/sigma));
s11 = exp(-abs(Dis11/sigma));
s21 = exp(-abs(Dis21/sigma));
s22 = exp(-abs(Dis22/sigma));

for i=1:89
   [vtm ptm]= sort(s11(:,i),'descend');s11(ptm(numK+1:end),i)=0;;
   [vtm ptm]= sort(s12(:,i),'descend');s12(ptm(Knear+1:end),i)=0;
   [vtm ptm]= sort(s21(:,i),'descend');s21(ptm(Knear+1:end),i)=0;
   [vtm ptm]= sort(s22(:,i),'descend');s22(ptm(numK +1:end),i)=0;
end
Wij1 = [s11.*ptmp11 s12.*ptmp;s21.*ptmp' s22.*ptmp22];
for i=1:178   
    if sum(Wij1(:,i))~=0
   Wij1(:,i) =  Wij1(:,i)/sum(Wij1(:,i));     
    end
end
% % % % Yl = eye(89);
% % % % % 从Ａ视角－－－－> B
 yu = inv(eye(89) - Wij1(90:end,90:end))*Wij1(90:end,1:89);
 
 yu2 = inv(eye(89) - Wij1(1:89,1:89))*Wij1(1:89,90:end);
for i=1:89
    if sum(yu(:,i))~=0
    [vtm ptm]=sort(yu(:,i),'descend'); yu(ptm(numK+1:end),i)= 0; 
    yu(:,i) = yu(:,i)/sum(yu(:,i));end
    if sum(yu2(:,i))~=0
    [vtm ptm]=sort(yu2(:,i),'descend'); yu2(ptm(numK+1:end),i)= 0; 
    yu2(:,i) = yu2(:,i)/sum(yu2(:,i));end
end
% % % % 
% % % % 
  POSS = yu.*yu2';
% % % % %   POSS = yu+yu2';

ppp = zeros(89,89); 
for i=1:89
    if max(POSS(:,i))~=0
  th= find(POSS(:,i)==max(POSS(:,i)));
  ppp(th,i)= 1;
    end
end

ppp2 = zeros(89,89); 
for i=1:89
    if max(POSS(i,:))~=0
  th= find(POSS(i,:)==max(POSS(i,:)));
   
  ppp2(th,i)= 1;
    end
end
ppp = ppp.*ppp2';
% % % %  ppp = ptmp;
fea_add_toU = [];y_add_toU = [];
for i=1:89
    if sum(ppp(:,i))==1
     
       fea_add_toU = [fea_add_toU fea_tr_a(:,find(tr_label_a==Uina(i))) fea_tr_b(:,find(tr_label_b==Uinb(find(ppp(:,i)==1))))]; 
        y_add_toU = [y_add_toU repmat(i,[1 length(find(tr_label_a==Uina(i)))+length(find(tr_label_b==Uinb(find(ppp(:,i)==1))))])];
    end
end
% % % % AddS(kk) = length(unique(y_add_toU));

% % % 
 [Wt3, Mt3,inCov3,exCov3] = XQDA_ICCV(fea_tr_a', fea_tr_b',tr_label_a,  tr_label_b, fea_add_toU',y_add_toU');
 M = Wt3*Mt3*Wt3';

  M1=M;M2=M;

 AddS(kk) = sum(ppp(:));

if kk>3
    if ((AddS(kk)==AddS(kk-1))  )|| kk>10
        break;
    end
end
kk=kk+1;

end

 
  clear AddS  
% % % tic
dist_M = MahDist(M, fea_te_b', fea_te_a');
for ii=1:89
    for jj=1:89    
    temp = dist_M(find(te_label_b==te_id(jj)),find(te_label_a==te_id(ii)));
    DistM(jj,ii)= sum(temp(:))/length(find(te_label_b==te_id(jj)))/length(find(te_label_a==te_id(ii)));
    end
end
toc

CMC(iter,:) =  EvalCMC(-DistM, 1:89,1:89,89);

end
% save('./final/PRID_result_gog_kis.mat','CMC','CMC_KIS','Var','Var_right')

an=0;cn=0;
for i=1:10
an = an + Var{i}(end);cn = cn + Var_right{i}(end);end

% CMC_K(iK,:) = mean(CMC);
 