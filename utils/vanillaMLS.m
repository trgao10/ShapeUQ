function [projPt] = vanillaMLS(pt, pCloud)
%VANILLAMLS Summary of this function goes here
%   Detailed explanation goes here

if (size(pt,1) == 1)
    pt = pt';
end

if (size(pCloud,1) > size(pCloud,2))
    pCloud = pCloud';
end

nV = size(pCloud,2);

%%% initial guess
distTestPt2PCloud = pdist2(pt',pCloud');
sd = sort(distTestPt2PCloud);
iSigma = sd(5);
Weights = exp(-distTestPt2PCloud.^2/iSigma^2);
coeff = pca(pCloud(:,Weights>1e-3)','Weights',Weights(Weights>1e-3));
EstNormal = coeff(:,3);

%%% non-linear optimization in one step

% Proj = @(r,t,normal) r+t*normal;
% sF = @(t) ((EstNormal'*(pCloud-repmat(Proj(pt,t,EstNormal),1,nV))).^2)*exp(-pdist2(Proj(pt,t,EstNormal)',pCloud').^2/iSigma^2)'/sum(exp(-pdist2(Proj(pt,t,EstNormal)',pCloud').^2/iSigma^2));
sF = @(t) ValueOnLine(pt, pCloud, t, EstNormal, iSigma);
t = fminbnd(sF,-iSigma,iSigma);
initProjPt = Proj(pt,t,EstNormal);

A = sort(squareform(pdist(pCloud'))+diag(Inf(nV,1)),2);
sigma = mean(A(:,10)); %% use the mean distance of each vertex to its 5-th neighbor as bandwidth
F = @(p) ValueAtPoint(pt, pCloud, p, sigma);
% options = optimoptions(@fminunc, 'Algorithm', 'trust-region', 'GradObj', 'on', 'Display', 'none');
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'GradObj', 'on', 'Display', 'none');
% options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'Display', 'none');
projPt = fminunc(F,initProjPt,options);
% F = @(p) ((pt-p)'*(pCloud-repmat(p,1,nV))).^2./((pt-p)'*(pt-p))*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2)');

end

function vLine = ValueOnLine(pt, pCloud, t, EstNormal, iSigma)
    nV = max(size(pCloud));
    vLine = ((EstNormal'*(pCloud-repmat(Proj(pt,t,EstNormal),1,nV))).^2)*SumNor(exp(-pdist2(Proj(pt,t,EstNormal)',pCloud').^2/iSigma^2))';%/sum(exp(-pdist2(Proj(pt,t,EstNormal)',pCloud').^2/iSigma^2));
end

% function vPt = ValueAtPoint(pt, pCloud, p, sigma)
%     nV = max(size(pCloud));
%     vPt = ((pt-p)'*(pCloud-repmat(p,1,nV))).^2./((pt-p)'*(pt-p))*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2)');
% end

function [vPt, gPt] = ValueAtPoint(pt, pCloud, p, sigma)
    nV = max(size(pCloud));
    vPt = ((pt-p)'*(pCloud-repmat(p,1,nV))).^2./((pt-p)'*(pt-p))*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2)');
    gPt = ((-2/sigma^2)*((repmat(p,1,nV)-pCloud-repmat((repmat(p,1,nV)-pCloud)*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2))',1,nV)).*repmat(((pt-p)'*(pCloud-repmat(p,1,nV))).^2./((pt-p)'*(pt-p)),3,1))...
        +2*(repmat((((pt-p)'*(pCloud-repmat(p,1,nV)))./((pt-p)'*(pt-p)).^2),3,1).*((pt-p)'*(pt-p)*(2*repmat(p,1,nV)-repmat(pt,1,nV)-pCloud)-repmat((pt-p)'*(pCloud-repmat(p,1,nV)),3,1).*repmat(p-pt,1,nV))))...
        *SumNor(exp(-pdist2(p',pCloud').^2/sigma^2)');
%     gPt = zeros(3,1);
%     for j=1:3
%         gPt(j) = ((-2/sigma^2)*((p(j)-pCloud(j,:)-(p(j)-pCloud(j,:))*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2))').*((pt-p)'*(pCloud-repmat(p,1,nV))).^2./((pt-p)'*(pt-p)))+2*((((pt-p)'*(pCloud-repmat(p,1,nV)))./((pt-p)'*(pt-p)).^2).*((pt-p)'*(pt-p)*(2*p(j)-pt(j)-pCloud(j,:))-(pt-p)'*(pCloud-repmat(p,1,nV))*(p(j)-pt(j)))))*SumNor(exp(-pdist2(p',pCloud').^2/sigma^2)');
%     end
end

function snv = SumNor(v)
    snv = v/sum(v);
end

function projPt = Proj(r, t, normal)
    projPt = r+t*normal;
end
