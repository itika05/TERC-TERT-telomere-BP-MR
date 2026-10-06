import numpy as np
from scipy.stats import norm
def pool(b,s):
    w=1/s**2; bf=np.sum(w*b)/w.sum(); sf=np.sqrt(1/w.sum()); Q=np.sum(w*(b-bf)**2); k=len(b)
    tau2=max(0,(Q-(k-1))/(w.sum()-np.sum(w**2)/w.sum())) if k>1 else 0
    wr=1/(s**2+tau2); br=np.sum(wr*b)/wr.sum(); sr=np.sqrt(1/wr.sum())
    from scipy.stats import chi2
    return dict(k=k,b_fixed=bf,se_fixed=sf,p_fixed=2*norm.sf(abs(bf/sf)),b_random=br,se_random=sr,p_random=2*norm.sf(abs(br/sr)),
                Q=Q,Q_p=chi2.sf(Q,k-1) if k>1 else np.nan,I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0,tau2=tau2)
