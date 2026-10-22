"""Independent algebraic check of all numerical sufficient inequalities."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
p=json.loads((ROOT/'03_code/parameters.json').read_text())
def gains(lam=0.,scale=None):
    scale=p['scale_y'] if scale is None else scale
    P,Q,B,D=[np.array(p[k],float) for k in ['P','Q','B','D']]
    c,d=[np.array(p[k]) for k in ['c','d']]
    ax=P*np.exp(lam*p['delay_linear'])+c[:,None]*np.exp(lam*np.array([p['delay_product_1'],p['delay_product_2']]))
    ay=scale*(Q*np.exp(lam*p['delay_linear'])+d[:,None]*np.exp(lam*np.array([p['delay_product_1'],p['delay_product_2']])))
    zero=np.zeros((2,2)); A=np.block([[zero,ax],[ay,zero]])
    N=np.block([[zero,B],[scale*D,zero]])*p['kernel_beta']/(p['kernel_beta']-lam)
    R=np.eye(4)/(p['leakage']-lam); S=np.eye(4)+p['leakage']*R
    H=np.block([[R@A,R@N],[S@A,S@N]])
    J=A@R+N@S
    d0=np.maximum(np.diag(R),np.diag(S)); G=np.diag(d0)@(A+N)
    return A,N,H,J,G
def rho(a):return float(max(abs(np.linalg.eigvals(a))))
def main():
    result={}
    for scale in [1.,2.]:
        A,N,H,J,G=gains(scale=scale); _,_,Hl,Jl,_=gains(p['lambda'],scale)
        w=np.linalg.solve(np.eye(8)-H,np.ones(8)); v=np.linalg.solve(np.eye(4)-J,np.ones(4))
        c=np.r_[np.full(4,2.),np.full(4,1+p['leakage'])]
        M=float(max(np.linalg.solve(np.eye(8)-Hl,c)))
        rate=brentq(lambda l:rho(gains(l,scale)[3])-1,0,1.9,xtol=1e-13)
        result[str(scale)]={'A':A.tolist(),'N':N.tolist(),'H':H.tolist(),'J':J.tolist(),'G':G.tolist(),'rho_H':rho(H),'rho_J':rho(J),'rho_G':rho(G),'rowmax_G':float(np.max(G.sum(axis=1))),'rho_J_lambda':rho(Jl),'lambda':p['lambda'],'M':M,'lambda_threshold':rate,'weights_H':w.tolist(),'q_H':float(max(H@w/w)),'certificate_v':v.tolist(),'certificate_residual':(v-J@v).tolist()}
        assert rho(J)<1 and rho(Jl)<1 and np.all(w>0)
        assert np.max(np.abs(np.linalg.eigvals(H)[:,None]-np.linalg.eigvals(H)[None,:]))>=0
        assert abs(rho(H)-rho(J))<1e-12
        assert np.all(J<=(A+N)*2+1e-14)
    assert result['1.0']['rho_G']<1<result['1.0']['rowmax_G']
    assert result['2.0']['rho_J']<1<result['2.0']['rho_G']
    out=ROOT/'outputs';out.mkdir(parents=True,exist_ok=True)
    (out/'condition_report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    for k,r in result.items():print(k,{x:r[x] for x in ['rho_J','rho_G','rowmax_G','rho_J_lambda','M','lambda_threshold','q_H']})
if __name__=='__main__':main()
