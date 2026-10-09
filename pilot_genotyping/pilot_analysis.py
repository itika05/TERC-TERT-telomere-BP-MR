"""Pakistani pilot genotyping: TERT rs2736100 and TERC rs10936599.
Genotype counts as reported in the July 2026 draft of the hTERT manuscript (TERC controls reconstructed
from the reported allelic OR 3.05 [1.106-8.445] and dominant OR 5.74; to be verified against the lab records).
Reference frequencies: Ensembl (1000 Genomes phase 3 PJL, SAS; gnomAD genomes SAS)."""
import json, numpy as np
from math import lgamma, exp, log, sqrt
from scipy import stats

def hwe_exact(nAA,nAB,nBB):
    n=nAA+nAB+nBB; nA=2*nAA+nAB; rare=min(nA,2*n-nA)
    mid=int(rare*(2*n-rare)/(2*n)); mid+= (mid%2)!=(rare%2)
    P={mid:1.0}; hr=(rare-mid)//2; hc=n-mid-hr
    h,a,b=mid,hr,hc
    while h>=2: P[h-2]=P[h]*h*(h-1)/(4*(a+1)*(b+1)); h-=2;a+=1;b+=1
    h,a,b=mid,hr,hc
    while h<=rare-2: P[h+2]=P[h]*4*a*b/((h+2)*(h+1)); h+=2;a-=1;b-=1
    t=sum(P.values()); obs=P[nAB]/t
    return min(1.0,sum(v/t for v in P.values() if v/t<=obs+1e-12))
def freeman_halton(tab):
    tab=np.array(tab); r=tab.sum(1); c=tab.sum(0); N=tab.sum()
    lp=lambda t: sum(lgamma(x+1) for x in r)+sum(lgamma(x+1) for x in c)-lgamma(N+1)-sum(lgamma(x+1) for x in t.flatten())
    o=lp(tab); p=0
    for x0 in range(min(r[0],c[0])+1):
        for x1 in range(min(r[0]-x0,c[1])+1):
            x2=r[0]-x0-x1
            if x2<0 or x2>c[2]: continue
            t=np.array([[x0,x1,x2],[c[0]-x0,c[1]-x1,c[2]-x2]])
            if (t<0).any(): continue
            l=lp(t)
            if l<=o+1e-9: p+=exp(l)
    return min(1,p)
def trend_test(case,ctrl):  # Cochran-Armitage, scores 0,1,2 (copies of the effect allele)
    s=np.array([0,1,2]); R=np.array(case); S=np.array(ctrl); N=R+S; n1=R.sum(); n=N.sum()
    T=(s*(R*S.sum()-S*n1)).sum()
    V=n1*S.sum()/n*(n*(s**2*N).sum()-((s*N).sum())**2)
    z=T/sqrt(V); return z, 2*stats.norm.sf(abs(z))
def firth_additive(case,ctrl):
    x=np.r_[np.repeat([0,1,2],case),np.repeat([0,1,2],ctrl)].astype(float); y=np.r_[np.ones(sum(case)),np.zeros(sum(ctrl))]
    X=np.c_[np.ones_like(x),x]
    def fit(fix=None):
        b=np.zeros(2); 
        if fix is not None: b[1]=fix
        for _ in range(200):
            p=1/(1+np.exp(-X@b)); W=p*(1-p); I=(X.T*W)@X; Ii=np.linalg.inv(I)
            h=np.einsum('ij,jk,ik->i',X*np.sqrt(W)[:,None],Ii,X*np.sqrt(W)[:,None])
            U=X.T@(y-p+h*(0.5-p))
            if fix is None: st=Ii@U
            else: st=np.array([U[0]/I[0,0],0.0])
            b=b+st
            if np.abs(st).max()<1e-10: break
        p=1/(1+np.exp(-X@b)); W=p*(1-p); I=(X.T*W)@X
        return b, np.sum(y*np.log(p)+(1-y)*np.log(1-p))+0.5*np.log(np.linalg.det(I))
    b,l=fit(); 
    def dev(v): return 2*(l-fit(v)[1])
    lo=b[1]; 
    while dev(lo)<3.8415: lo-=0.01
    hi=b[1]
    while dev(hi)<3.8415: hi+=0.01
    from scipy.optimize import brentq
    lo=brentq(lambda v:dev(v)-3.8415,lo,b[1]); hi=brentq(lambda v:dev(v)-3.8415,b[1],hi)
    return np.exp(b[1]), np.exp(lo), np.exp(hi), stats.chi2.sf(dev(0.0),1)
def power_allelic(maf,OR,ncase,nctrl,alpha=0.05,sims=4000,seed=1):
    rng=np.random.default_rng(seed); q0=maf; odds=q0/(1-q0)*OR; q1=odds/(1+odds); hit=0
    for _ in range(sims):
        a=rng.binomial(2*ncase,q1); c=rng.binomial(2*nctrl,q0)
        if stats.fisher_exact([[a,2*ncase-a],[c,2*nctrl-c]])[1]<alpha: hit+=1
    return hit/sims

ref={'rs10936599':{'allele':'T','PJL':(52,192),'SAS':(233,978),'gnomAD_SAS':(1299,4828)},
     'rs2736100':{'allele':'A','PJL':(71,192),'SAS':(386,978),'gnomAD_SAS':(1972,4800)}}
public={'rs10936599':{'LTL_EU':-0.1056,'HTN_bottomline_logOR':-0.0211,'HTN_bottomline_P':7.3e-7,'SBP_SA':(0.0015,0.89)},
        'rs2736100':{'LTL_EU':-0.0698,'HTN_EU_logOR':-0.0023,'HTN_EU_P':0.42,'SBP_SA':(-0.0096,0.33)}}
# genotype order: homozygous other / heterozygous / homozygous effect (T for TERC, A for TERT = LTL-shortening allele)
data={'rs10936599':{'gene':'TERC','case':[32,7,6],'ctrl':[3,5,2],'ctrl_source':'reconstructed from reported allelic and dominant ORs'},
      'rs2736100':{'gene':'TERT','case':[18,21,11],'ctrl':[5,2,3],'ctrl_source':'reported'}}
out={}
for rs,d in data.items():
    a=ref[rs]['allele']; r={}
    for g,cnt in (('cases',d['case']),('controls',d['ctrl']),('all',[x+y for x,y in zip(d['case'],d['ctrl'])])):
        n=sum(cnt); k=2*cnt[2]+cnt[1]; q=k/(2*n)
        pj=ref[rs]['PJL']
        r[g]=dict(n=n,genotypes=cnt,effect_allele=a,freq=round(q,3),het_obs=round(cnt[1]/n,3),het_exp=round(2*q*(1-q),3),
                  hwe_p=round(hwe_exact(*cnt),4),
                  vs_PJL_fisher_p=round(stats.fisher_exact([[k,2*n-k],[pj[0],pj[1]-pj[0]]])[1],3))
    ca,cn=2*d['case'][2]+d['case'][1],2*sum(d['case']); ka,kn=2*d['ctrl'][2]+d['ctrl'][1],2*sum(d['ctrl'])
    OR=(ca*(kn-ka))/((cn-ca)*ka); se=sqrt(1/ca+1/(cn-ca)+1/ka+1/(kn-ka))
    r['allelic_OR_effect_allele']=[round(OR,2),round(exp(log(OR)-1.96*se),2),round(exp(log(OR)+1.96*se),2)]
    r['allelic_fisher_p']=round(stats.fisher_exact([[ca,cn-ca],[ka,kn-ka]])[1],3)
    r['genotypic_freeman_halton_p']=round(freeman_halton([d['case'],d['ctrl']]),3)
    z,p=trend_test(d['case'],d['ctrl']); r['trend_p']=round(p,3)
    f=firth_additive(d['case'],d['ctrl']); r['firth_additive_OR']=[round(x,2) for x in f[:3]]+[round(f[3],3)]
    r['reference_freq']={k:round(v[0]/v[1],3) for k,v in ref[rs].items() if k!='allele'}
    r['public']=public[rs]; r['control_counts']=d['ctrl_source']
    r['power_allelic_OR1.5_2_3']=[power_allelic(ref[rs]['PJL'][0]/ref[rs]['PJL'][1],o,sum(d['case']),sum(d['ctrl'])) for o in (1.5,2,3)]
    out[rs]=r
json.dump(out,open('pilot_results.json','w'),indent=1)
print(json.dumps(out,indent=1))
