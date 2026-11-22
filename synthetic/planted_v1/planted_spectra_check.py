"""Planted-spectrum recovery tests for PIC Foundations (Objective 1 metrics + Objective 2 CUR rules).
SYNTHETIC ONLY. Matrices have KNOWN singular values; checks that each metric/rule recovers them.
Seed 6310. Shapes match the preprocessed synthetic_v1 train matrices (9505 x 82, 9505 x 59)."""
import numpy as np, scipy.linalg as sla, json
rng = np.random.default_rng(6310)

def planted(n, p, sig, noise=0.0):
    U,_ = np.linalg.qr(rng.standard_normal((n,p))); V,_ = np.linalg.qr(rng.standard_normal((p,p)))
    A = (U*sig) @ V.T
    if noise: A = A + noise*rng.standard_normal((n,p))
    return A

def metrics(A):
    s = np.linalg.svd(A, compute_uv=False)
    sr = (s**2).sum()/s[0]**2
    q = s/s.sum(); er = np.exp(-(q[q>0]*np.log(q[q>0])).sum())
    k95_2 = int(np.searchsorted(np.cumsum(s**2)/(s**2).sum(), .95)+1)
    k95_1 = int(np.searchsorted(np.cumsum(s)/s.sum(), .95)+1)
    i = np.arange(1, len(s)+1); m = slice(0, min(30, len(s)))
    pl = np.polyfit(np.log(i[m]), np.log(s[m]), 1); ex = np.polyfit(i[m], np.log(s[m]), 1)
    rpl = np.sum((np.log(s[m])-np.polyval(pl, np.log(i[m])))**2); rex = np.sum((np.log(s[m])-np.polyval(ex, i[m]))**2)
    # permutation null on column-permuted copy
    P = np.column_stack([rng.permutation(c) for c in A.T]); sp = np.linalg.svd(P, compute_uv=False)
    n,p = A.shape; Z=(A-A.mean(0))/A.std(0); sz=np.linalg.svd(Z,compute_uv=False)/np.sqrt(n)
    mp = 1+np.sqrt(p/n)
    return dict(stable_rank=round(sr,2), entropy_erank=round(er,2), k95_sigma2=k95_2, k95_sigma=k95_1,
                decay_best="power" if rpl<rex else "exponential", alpha=round(-pl[0],3), beta=round(-ex[0],3),
                above_perm_null=int((s>sp[0]).sum()), above_MP_edge=int((sz>mp).sum()), kappa=round(s[0]/s[-1],1))

def cur_ratio(A, k, rule):
    U,s,Vt = np.linalg.svd(A, full_matrices=False); err_opt = np.sqrt((s[k:]**2).sum())
    if rule=="cpqr":
        J = sla.qr(A, pivoting=True, mode='r')[1][:k]; I = sla.qr(A.T, pivoting=True, mode='r')[1][:k]
    elif rule=="deim":
        def deim(W):
            idx=[int(np.argmax(abs(W[:,0])))]
            for j in range(1,k):
                c=np.linalg.solve(W[idx,:j], W[idx,j]); r=W[:,j]-W[:,:j]@c; idx.append(int(np.argmax(abs(r))))
            return np.array(idx)
        J=deim(Vt[:k].T); I=deim(U[:,:k])
    else:
        lv=(Vt[:k]**2).sum(0); lu=(U[:,:k]**2).sum(1)
        if rule=="lev_det": J=np.argsort(-lv)[:k]; I=np.argsort(-lu)[:k]
        elif rule=="lev_rand": J=rng.choice(A.shape[1],k,False,lv/lv.sum()); I=rng.choice(A.shape[0],k,False,lu/lu.sum())
        else: J=rng.choice(A.shape[1],k,False); I=rng.choice(A.shape[0],k,False)
    C=A[:,J]; R=A[I,:]; Um=np.linalg.pinv(C)@A@np.linalg.pinv(R)
    return np.linalg.norm(A-C@Um@R)/err_opt

out={}
for name,(n,p) in {"Xeng_shape":(9505,82),"Xraw_shape":(9505,59)}.items():
    i=np.arange(1,p+1)
    cases={"exact_rank10_plus_noise":(np.r_[np.linspace(20,10,10),np.zeros(p-10)],0.01),
           "power_alpha1.5":(10*i**-1.5,0),"exponential_beta0.15":(10*np.exp(-0.15*i),0)}
    for cname,(sig,noise) in cases.items():
        A=planted(n,p,sig,noise); r=metrics(A)
        k=10; r["CUR_ratio_k10_median"]={rule: round(float(np.median([cur_ratio(A,k,rule) for _ in range(5 if rule in('uniform','lev_rand') else 1)])),3)
                                         for rule in ["uniform","lev_rand","lev_det","cpqr","deim"]}
        out[f"{name}/{cname}"]=r; print(f"{name}/{cname}", json.dumps(r))
json.dump(out, open("planted_results.json","w"), indent=1)
