function feat = gen_10_iter(feature,numPCA,feature_path)
view_a_feature =  importdata('/media/lzm/My Passport/lzm/view_a_feature_ILIDS.txt');
view_b_feature =  importdata('/media/lzm/My Passport/lzm/view_b_feature_ILIDS.txt');
view_a_id =  importdata('/media/lzm/My Passport/lzm/view_a_id_ILIDS.txt');
view_b_id =  importdata('/media/lzm/My Passport/lzm/view_b_id_ILIDS.txt');
% view_a_feature  = feature.view_a_feature;
% view_b_feature  = feature.view_b_feature;
% view_a_id  = feature.view_a_id;
% view_b_id  = feature.view_b_id;
% half train, half test {OK for viper, ilids, prid2011, prid450s, cuhk01}
view_a_feature=view_a_feature';
view_b_feature=view_b_feature';
view_a_id=view_a_id';
view_b_id=view_b_id';
total_person = length(unique(view_a_id));
dim_fea = size(view_a_feature,1);
% dimPCA = 100:100:1000;
% for numP = 1:10
numPCA = 600;
rand('state',1)
% Image_SET = importdata(['/media/lzm/000D077800061AC2/Re-ID/DLcode/semi-transfer/stage1/label/PRID_halftrain_iter' num2str(iter_num) '_imgs_truelabel.txt']);
% total_img_num =  length(Image_SET.textdata);
%     numPCA = dimPCA(numP);
for iter = 1:10
  %% partition the feature into train and test
     p = randperm(total_person);
     tr_id = p(1:floor(total_person/2));
     te_id = p(1+floor(total_person/2):end);
     tr_label_a=[];tr_label_b=[];te_label_a=[];te_label_b=[];
     % divide label
     for ii=1:length(tr_id)
      tr_label_a = [tr_label_a  view_a_id(find(view_a_id == tr_id(ii)))];
      tr_label_b =  [tr_label_b  view_b_id(find(view_b_id == tr_id(ii)))];
     end
     for ii=1:length(te_id)
      te_label_a = [te_label_a  view_a_id(find(view_a_id == te_id(ii)))];
      te_label_b = [te_label_b  view_b_id(find(view_b_id == te_id(ii)))];
     end
     % construct feature space for train & test
     fea_tr_a = zeros(dim_fea,length(tr_label_a));
     fea_tr_b = zeros(dim_fea,length(tr_label_b));
     fea_te_a = zeros(dim_fea,length(te_label_a));
     fea_te_b = zeros(dim_fea,length(te_label_b));    
     % divide feature
     Luni_a = tr_id;
     Luni_b = tr_id;
     Luni_ae = te_id;
     Luni_be =te_id;

     for ii=1:150
     fea_tr_a(:,find(Luni_a(ii)==tr_label_a)) =  view_a_feature(:,find(view_a_id == tr_id(ii)));
     fea_tr_b(:,find(Luni_b(ii)==tr_label_b)) =  view_b_feature(:,find(view_b_id == tr_id(ii)));
     fea_te_a(:,find(Luni_ae(ii)==te_label_a)) =  view_a_feature(:,find(view_a_id == te_id(ii)));
     fea_te_b(:,find(Luni_be(ii)==te_label_b)) =  view_b_feature(:,find(view_b_id == te_id(ii)));
     end
  %% PCA reduce dim to numPCA
     X = [fea_tr_a fea_tr_b] ;
     [W, S] = PCA(X', numPCA);
     
     m = mean(X,2);
     fea_tr_a= S(1:size(fea_tr_a,2),1:numPCA)';
     fea_tr_b= S(1+size(fea_tr_a,2):end,1:numPCA)';
     fea_te_a= W(:,1:numPCA)'*fea_te_a;
     fea_te_b= W(:,1:numPCA)'*fea_te_b;
%      fea_tr_a=u(:,1:numPCA)'*(fea_tr_a-repmat(m,1,size(fea_tr_a,2)));
%      fea_tr_b=u(:,1:numPCA)'*(fea_tr_b-repmat(m,1,size(fea_tr_b,2)));
%      fea_te_a=u(:,1:numPCA)'*(fea_te_a-repmat(m,1,size(fea_te_a,2)));
%      fea_te_b=u(:,1:numPCA)'*(fea_te_b-repmat(m,1,size(fea_te_b,2)));  
   %% save feature
     feat{iter}.fea_tr_a = fea_tr_a;
     feat{iter}.fea_tr_b = fea_tr_b;
     feat{iter}.fea_te_a = fea_te_a;
     feat{iter}.fea_te_b = fea_te_b;
     feat{iter}.tr_label_a = tr_label_a;
     feat{iter}.tr_label_b = tr_label_b;
     feat{iter}.te_label_a = te_label_a;
     feat{iter}.te_label_b = te_label_b;
     feat{iter}.tr_id = tr_id;
     feat{iter}.te_id = te_id;     
end
%  save([feature_path '_10_iteration.mat'],'feat','-v7.3');
save(['ILIDS_GOGdim' num2str(numPCA) '_10_iteration.mat'],'feat','-v7.3');
% numP
end
% end

