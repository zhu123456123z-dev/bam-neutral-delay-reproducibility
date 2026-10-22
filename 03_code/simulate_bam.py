"""Reproducible method-of-steps solver, independent of the MATLAB driver.
Each integration interval is at most the minimum discrete delay. Delayed
values come only from already completed dense-output segments.
"""
from pathlib import Path
import os,json,bisect,time
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mpl'))
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.io import savemat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=json.loads((ROOT/'03_code/parameters.json').read_text());out=ROOT/'outputs'/'python'
P,Q,B,D=[np.array(p[k]) for k in ['P','Q','B','D']]
c,d=[np.array(p[k]) for k in ['c','d']];a,b,e=[np.array(p[k]) for k in ['input_sin','input_cos','input_ergodic']]
lags=[p[k] for k in ['delay_linear','delay_product_1','delay_product_2']]

def solve(t0,t1,history,fine=False):
    bounds=[float(t0)];segments=[];state=np.r_[history,np.zeros(4)]
    def delayed(t):
        if t<=t0+1e-12:return np.r_[history,np.zeros(4)]
        ix=min(bisect.bisect_left(bounds,t)-1,len(segments)-1)
        if ix<0:raise RuntimeError('Uncomputed history')
        assert t<=bounds[-1]+2e-10,(t,bounds[-1])
        return segments[ix](min(t,bounds[-1]))
    def rhs(t,z):
        Z=[delayed(t-delay) for delay in lags]
        u=a*np.sin(t)+b*np.cos(np.sqrt(2)*t)+e/(1+t*t)
        dx=-p['leakage']*z[:2]+(.9+.1*np.sin(t))*P@np.tanh(Z[0][2:4])+B@z[6:8]+c*np.tanh(Z[1][2])*np.tanh(Z[2][3])+u[:2]
        dy=-p['leakage']*z[2:4]+p['scale_y']*((.9+.1*np.cos(np.sqrt(2)*t))*Q@np.tanh(Z[0][:2])+D@z[4:6]+d*np.tanh(Z[1][0])*np.tanh(Z[2][1]))+u[2:4]
        dz=np.r_[dx,dy];return np.r_[dz,p['kernel_beta']*(np.tanh(dz)-z[4:])]
    steps=int(round((t1-t0)/min(lags))); factor=10 if fine else 1
    for k in range(steps):
        left=t0+k*min(lags);right=min(t0+(k+1)*min(lags),t1)
        sol=solve_ivp(rhs,(left,right),state,method='DOP853',rtol=p['relative_tolerance']/factor,atol=p['absolute_tolerance']/factor,max_step=p['max_step']/(2 if fine else 1),dense_output=True)
        assert sol.success,sol.message
        segments.append(sol.sol);bounds.append(right);state=sol.y[:,-1]
    def evaluate(ts, derivative=False):
        values=np.column_stack([segments[max(0,min(bisect.bisect_left(bounds,float(t))-1,len(segments)-1))](float(t)) for t in ts])
        if derivative:
            # Evaluate the model RHS; at t0 this is the right-hand derivative.
            return np.column_stack([rhs(float(t),values[:,j])[:4] for j,t in enumerate(ts)])
        return values
    return evaluate

def plot_derivative_errors(t, derivative_errors, output, decay_rate):
    """Plot E_D(t), with an exponential slope guide (not a theorem bound)."""
    fig,axs=plt.subplots(1,2,figsize=(10.5,4))
    for k,error in enumerate(derivative_errors):
        axs[0].plot(t,error,label=f'history {k+1}')
        axs[1].semilogy(t,np.maximum(error,1e-16),label=f'history {k+1}')
    guide=np.max(derivative_errors[:,0])*np.exp(-decay_rate*t)
    axs[1].semilogy(t,guide,'k--',label=rf'slope guide $e^{{-{decay_rate:g}t}}$')
    axs[0].set(xlim=(0,min(15,t[-1])),xlabel='t',ylabel='Derivative sup error')
    axs[1].set(xlim=(0,min(25,t[-1])),ylim=(1e-15,max(10,2*np.max(derivative_errors))),xlabel='t',ylabel='Derivative sup error')
    for ax in axs:ax.grid(alpha=.2)
    axs[1].legend(fontsize=9);fig.tight_layout()
    fig.savefig(output/'fig5_derivative_errors.png',dpi=300)
    fig.savefig(output/'fig5_derivative_errors.pdf')
    plt.close(fig)

def main():
    out.mkdir(parents=True,exist_ok=True)
    t=np.arange(0,p['end_time']+p['output_step']/2,p['output_step']);start=time.time()
    reference=solve(-p['burn_in'],p['end_time'],np.zeros(4))
    zr=reference(t);reference_derivatives=reference(t,derivative=True);print('reference done',flush=True)
    refined=solve(-p['burn_in'],p['end_time'],np.zeros(4),True)
    zf=refined(t);zfp=refined(t,derivative=True);print('fine done',flush=True)
    longer=solve(-p['burn_in_check'],p['end_time'],np.zeros(4),True)
    zl=longer(t);zlp=longer(t,derivative=True);print('long done',flush=True)
    histories=np.array(p['histories']);trajectories=[];trajectory_derivatives=[]
    for h in histories:
        response=solve(0,p['end_time'],h)
        trajectories.append(response(t)[:4]);trajectory_derivatives.append(response(t,derivative=True))
    trajectories=np.array(trajectories);trajectory_derivatives=np.array(trajectory_derivatives)
    fine_response=solve(0,p['end_time'],histories[0],True)
    tf=fine_response(t)[:4];tfp=fine_response(t,derivative=True)
    errors=np.max(np.abs(trajectories-zr[:4]),axis=1)
    derivative_errors=np.max(np.abs(trajectory_derivatives-reference_derivatives),axis=1)
    joint_errors=np.maximum(errors,derivative_errors)
    validation={'solver':'Python SciPy DOP853 method of steps','scipy_version':scipy.__version__,'reference_refinement':float(np.max(abs(zr[:4]-zf[:4]))),'burn_in_difference':float(np.max(abs(zf[:4]-zl[:4]))),'trajectory_refinement':float(np.max(abs(trajectories[0]-tf))),'error_t10':errors[:,500].tolist(),'error_t20':errors[:,1000].tolist(),'elapsed_seconds':time.time()-start}
    validation.update(reference_derivative_refinement=float(np.max(abs(reference_derivatives-zfp))),
                      derivative_burn_in_difference=float(np.max(abs(zfp-zlp))),
                      trajectory_derivative_refinement=float(np.max(abs(trajectory_derivatives[0]-tfp))),
                      derivative_initial_errors=derivative_errors[:,0].tolist(),
                      derivative_error_t10=derivative_errors[:,500].tolist(),
                      derivative_error_t20=derivative_errors[:,1000].tolist())
    (out/'python_validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    np.savez(out/'simulation_python.npz',t=t,zr=zr,trajectories=trajectories,errors=errors,
             reference_derivatives=reference_derivatives,trajectory_derivatives=trajectory_derivatives,
             derivative_errors=derivative_errors,joint_errors=joint_errors)
    np.savetxt(out/'trajectories_python.csv',np.c_[t,zr[:4].T,errors.T],delimiter=',',header='t,x1,x2,y1,y2,error1,error2,error3',comments='')
    derivative_header='t,dx1_ref,dx2_ref,dy1_ref,dy2_ref,'+','.join(f'derivative_error{k+1}' for k in range(len(histories)))
    np.savetxt(out/'derivative_errors_python.csv',np.c_[t,reference_derivatives.T,derivative_errors.T],delimiter=',',header=derivative_header,comments='')
    names=[r'$x_1$',r'$x_2$',r'$y_1$',r'$y_2$']
    fig,axs=plt.subplots(2,2,figsize=(10.5,6.2))
    for j,ax in enumerate(axs.flat):ax.plot(t,zr[j],lw=1.2);ax.set(xlabel='t',ylabel=names[j]);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'fig1_states.png',dpi=220);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(10.5,6.2))
    for j,ax in enumerate(axs.flat):
        for k in range(3):ax.plot(t,trajectories[k,j],lw=1,label=f'history {k+1}')
        ax.plot(t,zr[j],'k--',lw=1.3,label='reference');ax.set(xlim=(0,15),xlabel='t',ylabel=names[j]);ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=9);fig.tight_layout();fig.savefig(out/'fig2_histories.png',dpi=220);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4))
    for k in range(3):
        axs[0].plot(t,errors[k],label=f'history {k+1}')
        axs[1].semilogy(t,np.maximum(errors[k],1e-16),label=f'history {k+1}')
    axs[1].semilogy(t,max(errors[:,0])*np.exp(-p['lambda']*t),'k--',label=r'slope guide $e^{-0.2t}$')
    axs[0].set(xlim=(0,15),xlabel='t',ylabel='State sup error');axs[1].set(xlim=(0,25),ylim=(1e-15,10),xlabel='t',ylabel='State sup error')
    for ax in axs:ax.grid(alpha=.2)
    axs[1].legend(fontsize=9);fig.tight_layout();fig.savefig(out/'fig3_errors.png',dpi=220);plt.close(fig)
    plot_derivative_errors(t,derivative_errors,out,p['lambda'])
    print(json.dumps(validation,indent=2),flush=True)
if __name__=='__main__':main()
