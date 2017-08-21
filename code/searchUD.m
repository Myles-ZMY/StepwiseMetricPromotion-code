function [fea_add_toU y_add_toU y_add_toU_real]= searchUD(fea_tr_a,fea_tr_b,tr_label_a,tr_label_b,M1,M2) 
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
LuniA = unique(tr_label_a,'stable');LuniB = unique(tr_label_b,'stable');
% dis11 = MahDist(M1, fea_tr_a', fea_tr_a');
% dis12 = MahDist(M1, fea_tr_a', fea_tr_b');
% dis21 = MahDist(M2, fea_tr_b', fea_tr_a');
% dis22 = MahDist(M2, fea_tr_b', fea_tr_b');
dis11 = MahDist(M1, fea_tr_a', fea_tr_a')+MahDist(M2, fea_tr_a', fea_tr_a') ;
dis21 =  MahDist(M1, fea_tr_b', fea_tr_a')+MahDist(M2, fea_tr_b', fea_tr_a') ;
dis12 =  MahDist(M2, fea_tr_a', fea_tr_b')+MahDist(M1, fea_tr_a', fea_tr_b') ;
dis22 = MahDist(M2, fea_tr_b', fea_tr_b')+ MahDist(M1, fea_tr_b', fea_tr_b');
numK = 10;

% % parpool open
% parfor ii=1:length(LuniA)
%     parfor jj=1:length(LuniB)
% %       dist_maxpooling2(jj,ii) = MahDist(M2, max(fea_te_b(:,find(te_label_b==te_id(jj))),[],2)'*W2,  max(fea_te_a(:,find(te_label_a==te_id(ii))),[],2)'*W2);
%         temp11 = dis11(find(tr_label_a==LuniA(jj)),find(tr_label_a==LuniA(ii)));
%         temp12 = dis12(find(tr_label_b==LuniB(jj)),find(tr_label_a==LuniA(ii)));
%         temp21 = dis21(find(tr_label_a==LuniA(ii)),find(tr_label_b==LuniB(jj)));
%         temp22 = dis22(find(tr_label_b==LuniB(ii)),find(tr_label_b==LuniB(jj)));
%         Dis11(jj,ii)= sum(temp11(:))/length(find(tr_label_a==LuniA(jj)))/length(find(tr_label_a==LuniA(ii)));
%         Dis12(jj,ii)= sum(temp12(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_a==LuniA(ii)));
%         Dis21(jj,ii)= sum(temp21(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_a==LuniA(ii)));
%         Dis22(jj,ii)= sum(temp22(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_b==LuniB(ii)));
%     end
% end
% parpool();
 for ii = 1:length(LuniA)
 Dis11(:,ii)= mean(reshape(mean(dis11(:,find(tr_label_a==LuniA(ii))),2),[20 length(LuniA)]),1);
 Dis21(:,ii)= mean(reshape(mean(dis21(:,find(tr_label_a==LuniA(ii))),2),[20 length(LuniB)]),1);
 end
 for ii = 1:length(LuniB)
 Dis12(:,ii)= mean(reshape(mean(dis12(:,find(tr_label_b==LuniB(ii))),2),[20 length(LuniA)]),1);
 Dis22(:,ii)= mean(reshape(mean(dis22(:,find(tr_label_b==LuniB(ii))),2),[20 length(LuniB)]),1);
 end
% parpool close
% poss1 = zeros(length(tr_label_b),length(tr_label_a));
% poss2 = zeros(length(tr_label_b),length(tr_label_a));
% YY_poss = zeros(length(tr_label_b),length(tr_label_a));
fea_addto_U = [];y_addto_U = [];
numK = 10;
% % % for ii=1:length(LuniA)
% % %     for jj=1:length(LuniB)
% % % %       dist_maxpooling2(jj,ii) = MahDist(M2, max(fea_te_b(:,find(te_label_b==te_id(jj))),[],2)'*W2,  max(fea_te_a(:,find(te_label_a==te_id(ii))),[],2)'*W2);        
% % %         temp21 = dis21(find(tr_label_b==LuniB(jj)),find(tr_label_a==LuniA(ii)));
% % %         temp12 = dis12(find(tr_label_a==LuniA(ii)),find(tr_label_b==LuniB(jj)));       
% % %         Dis12(jj,ii)= sum(temp12(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_a==LuniA(ii)));
% % %         Dis21(ii,jj)= sum(temp21(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_a==LuniA(ii)));
% % %         
% % %     end
% % % end
% % % 
% % % for ii=1:length(LuniA)
% % %     for jj=1:length(LuniA)
% % %          temp11 = dis11(find(tr_label_a==LuniA(jj)),find(tr_label_a==LuniA(ii)));
% % %          Dis11(jj,ii)= sum(temp11(:))/length(find(tr_label_a==LuniA(jj)))/length(find(tr_label_a==LuniA(ii)));
% % %     end
% % % end 
% % % 
% % % 
% % % 
% % % for ii=1:length(LuniB)
% % %     for jj=1:length(LuniB)
% % %         temp22 = dis22(find(tr_label_b==LuniB(ii)),find(tr_label_b==LuniB(jj)));
% % %         Dis22(jj,ii)= sum(temp22(:))/length(find(tr_label_b==LuniB(jj)))/length(find(tr_label_b==LuniB(ii)));
% % %  
% % %     end
% % % end 
for i=1:length(LuniA)
   [vtmp11(:,i) ptmp11(:,i)] = sort(Dis11(:,i));   
   [vtmp21(:,i) ptmp21(:,i)] = sort(Dis21(:,i));    
end
for i=1:length(LuniB)
   [vtmp12(:,i) ptmp12(:,i)] = sort(Dis12(:,i));
   [vtmp22(:,i) ptmp22(:,i)] = sort(Dis22(:,i)); 
end
s12 = exp(-abs(Dis12/3000));
s11 = exp(-abs(Dis11/3000));
s21 = exp(-abs(Dis21/3000));
s22 = exp(-abs(Dis22/3000));

for i=1:length(LuniA)
   [vtm ptm]= sort(s11(:,i),'descend');s11(ptm(numK+1:end),i)=0;   
   [vtm ptm]= sort(s21(:,i),'descend');s21(ptm(numK+1:end),i)=0;   
end
for i=1:length(LuniB)
   [vtm ptm]= sort(s12(:,i),'descend');s12(ptm(numK+1:end),i)=0;
   [vtm ptm]= sort(s22(:,i),'descend');s22(ptm(numK+1:end),i)=0;
end
Wij1 =   [s11 s12;s21 s22];
for i=1:length(LuniA)+length(LuniB)   
   Wij1(:,i) =  Wij1(:,i)/sum(Wij1(:,i));
end
Yl11 = bsxfun(@eq,LuniA,LuniA');Yl22=bsxfun(@eq,LuniB,LuniB');
yu = inv(eye(length(LuniB)) - Wij1(length(LuniA)+1:end,length(LuniA)+1:end))*Wij1(length(LuniA)+1:end,1:length(LuniA))*Yl11;
yu2 = inv(eye(length(LuniA)) - Wij1(1:length(LuniA),1:length(LuniA)))*Wij1(1:length(LuniA),length(LuniA)+1:end)*Yl22;
for i=1:length(LuniA)
    [vtm ptm]=sort(yu(:,i),'descend'); 
  yu(ptm(numK+1:end),i)= 0; 
    yu(:,i) = yu(:,i)/sum(yu(:,i));
end
for i=1:length(LuniB)
    [vtm ptm]=sort(yu2(:,i),'descend'); 
 yu2(ptm(numK+1:end),i)= 0; 
    yu2(:,i) = yu2(:,i)/sum(yu2(:,i));
end

POSS = yu+yu2';
Uina = unique(tr_label_a,'stable');
Uinb = unique(tr_label_b,'stable');
ppp = zeros(length(LuniB),length(LuniA)); 
vv=min(length(Uinb),length(Uina));

for i=1:length(LuniA)
 if max(POSS(:,i))~=0
  th= find(POSS(:,i)==max(POSS(:,i)));
  thrl = find(POSS(th,:)==max(POSS(th,:)));
  if thrl==i
  ppp(th,i)= 1;end
 end
end


ppp2 = zeros(length(LuniA),length(LuniB)); 
 for i=1:length(LuniB)
    if max(POSS(i,:))~=0
  th= find(POSS(i,:)==max(POSS(i,:)));
  thrl = find(POSS(:,th)==max(POSS(:,th)));
  if thrl==i
  ppp2(th,i)= 1;
  end
    end
 end
ppp  = ppp.*ppp2';

 if vv==length(Uina)

   
fea_add_toU = [];y_add_toU = [];y_add_toU_real=[];
for i=1:length(LuniA)
    if sum(ppp(:,i))==1
       fea_add_toU = [fea_add_toU fea_tr_a(:,find(tr_label_a==Uina(i))) fea_tr_b(:,find(tr_label_b==Uinb(find(ppp(:,i)==1))))]; 
        y_add_toU = [y_add_toU repmat(Uina(i),[1 length(find(tr_label_a==Uina(i)))+length(find(tr_label_b==Uinb(find(ppp(:,i)==1))))])];
     y_add_toU_real = [y_add_toU_real repmat(Uina(i),[1 length(find(tr_label_a==Uina(i)))])  repmat(Uinb(find(ppp(:,i)==1)),[1 length(find(tr_label_b==Uinb(find(ppp(:,i)==1))))])];
    end
end

else
fea_add_toU = [];y_add_toU = [];y_add_toU_real=[];
for i=1:length(LuniB)
    if sum(ppp(i,:))==1
       fea_add_toU = [fea_add_toU fea_tr_a(:,find(tr_label_a==Uina(find(ppp(i,:)==1)))) fea_tr_b(:,find(tr_label_b==Uinb(i)))]; 
        y_add_toU = [y_add_toU repmat(Uinb(i),[1 length(find(tr_label_a==Uina(find(ppp(i,:)==1))))+length(find(tr_label_b==Uinb(i)))])];
                 y_add_toU_real = [y_add_toU_real repmat(Uinb(i),[1 length(find(tr_label_b==Uinb(i)))])  repmat(Uina(find(ppp(i,:)==1)),[1 length(find(tr_label_a==Uina(find(ppp(i,:)==1))))])];
   
    end
end   
end
%  AddS(kk) = length(unique(y_add_toU));



% end
 

