function [W, M, inCov, exCov] = XQDA(galX, probX, galLabels, probLabels, V1X,V1L,options)
%% function [W, M, inCov, exCov] = XQDA(galX, probX, galLabels, probLabels, options)
% Cross-view Quadratic Discriminant Analysis for subspace and metric
% learning
%
% Input:
%   <galX>: features of gallery samples. Size: [n, d]
%   <probX>: features of probe samples. Size: [m, d]
%   <galLabels>: class labels of the gallery samples
%   <probLabels>: class labels of the probe samples
%   [options]: optional parameters. A structure containing any of the
%   following fields:
%       lambda: the regularizer. Default: 0.001
%       qdaDims: the number of dimensions to be preserved in the learned
%       subspace. Negative values indicate automatic dimension selection by
%       perserving latent values larger than 1. Default: -1.
%       verbose: whether to print the learning details. Default: false
%
% Output:
%   W: the subspace projection matrix. Size: [d, r], where r is the
%   subspace dimension.
%   M: the learned metric kernel. Size: [r,r]
%   inCov: covariance matrix of the intra-personal difference class. Size:
%   [r,r]
%   exCov: covriance matrix of the extra-personal difference class. Size:
%   [r,r]
%
%   With W, inCov, and exCov, we can quickly learn another metric of different dimensions:
%       W2 = W(:, 1:r2);
%       M2 = inv(inCov(1:r2, 1:r2)) - inv(exCov(1:r2, 1:r2));
% 
% Example:
%     Please see Demo_XQDA.m.
%
% Reference:
%   Shengcai Liao, Yang Hu, Xiangyu Zhu, and Stan Z. Li. Person
%   re-identification by local maximal occurrence representation and metric
%   learning. In IEEE Conference on Computer Vision and Pattern Recognition, 2015.
% 
% Version: 1.0
% Date: 2015-04-30
%
% Author: Shengcai Liao
% Institute: National Laboratory of Pattern Recognition,
%   Institute of Automation, Chinese Academy of Sciences
% Email: scliao@nlpr.ia.ac.cn
lambda = 0.001;
qdaDims = -1;
verbose = false;
 
[inCov exCov]=coutINEX(galX, galX, galLabels,galLabels);
%% v1
if length(unique(V1L))>1
[in2 ex2]=coutINEX(V1X,V1X,V1L,V1L);
else
in2 = 0;
ex2 = 0;
end
[in3 ex3]=coutINEX(probX, probX, probLabels,probLabels);
% [in3 en3]=coutINEX(V2X,V2X,V2L,V2L);
 inCov = inCov +  2*in2 +in3 ;
 exCov = exCov +  2*ex2 +ex3 ;

t0 = tic;
  [V, S] = svd(inCov \ exCov);

if verbose == true
    fprintf(' %.3g seconds.\n', toc(t0));
end

latent = diag(S);
[latent, index] = sort(latent, 'descend');
energy = sum(latent);
minv = latent(end);

r = sum(latent > 1);
energy = sum(latent(1:r)) / energy;

if qdaDims > r
    qdaDims = r;
end

if qdaDims <= 0
    qdaDims = max(1,r);
end

if verbose == true
    fprintf('Energy remained: %f, max: %f, min: %f, all min: %f, #opt-dim: %d, qda-dim: %d.\n', energy, latent(1), latent(max(1,r)), minv, r, qdaDims);
end

V = V(:, index(1:qdaDims));
 
if ~exist('W', 'var');
    W = V;
else
    W = W * V;
end

if verbose == true
    fprintf('Compute kernel matrix...');
end

t0 = tic;


  inCov = V' * inCov * V;
  exCov = V' * exCov * V;
  M = inv(inCov) - inv(exCov);
% % % 
% % % if verbose == true
% % %     fprintf(' %.3g seconds.\n\n', toc(t0));
% % % end
