import csv, glob, math, itertools, os, json
_repo = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
# folder with the result files of 05_simulated_analysis and 06_real_datasets: set EEGDN_RESULTS, or results/ inside or next to the repository
RD = os.environ.get('EEGDN_RESULTS') or next((d for d in (os.path.join(_repo, 'results'), os.path.join(os.path.dirname(_repo), 'results')) if os.path.isdir(d)), os.path.join(_repo, 'results'))
if not RD.endswith(os.sep): RD += os.sep
def gammainc_q(a,x):  # upper regularized gamma Q(a,x)
    if x<=0: return 1.0
    if x<a+1:
        s=1/a; t=s; n=a
        for _ in range(2000):
            n+=1; t*=x/n; s+=t
            if abs(t)<abs(s)*1e-15: break
        return 1-s*math.exp(-x+a*math.log(x)-math.lgamma(a))
    b=x+1-a; c=1e300; d=1/b; h=d
    for i in range(1,2000):
        an=-i*(i-a); b+=2; d=an*d+b; d=1e-300 if abs(d)<1e-300 else d
        c=b+an/c; c=1e-300 if abs(c)<1e-300 else c
        d=1/d; de=d*c; h*=de
        if abs(de-1)<1e-15: break
    return math.exp(-x+a*math.log(x)-math.lgamma(a))*h
def chi2_sf(x,k): return gammainc_q(k/2,x/2)
def ranks(v):
    o=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v); i=0; ties=[]
    while i<len(v):
        j=i
        while j+1<len(v) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        ties.append(j-i+1); i=j+1
    return r,ties
def friedman(M):  # M: n blocks x k treatments
    n=len(M); k=len(M[0]); R=[[0]*k for _ in range(n)]; T=0
    sums=[0]*k
    for i,row in enumerate(M):
        r,ties=ranks(row)
        for j in range(k): sums[j]+=r[j]
        T+=sum(t**3-t for t in ties)
    num=12*sum(s*s for s in sums)/(n*k*(k+1))-3*n*(k+1)
    den=1-T/(n*k*(k*k-1))
    if den<=0: return 1.0
    return chi2_sf(num/den,k-1)
def wilcoxon(x,y):
    d=[a-b for a,b in zip(x,y) if a!=b]; n=len(d)
    if n==0: return 1.0
    r,ties=ranks([abs(v) for v in d]); W=sum(ri for ri,v in zip(r,d) if v>0); tot=n*(n+1)/2
    if n<=20:
        cnt=0;N=2**n
        # exact via DP on doubled ranks
        rr=[int(round(2*v)) for v in r]; S=sum(rr); dp={0:1}
        for v in rr:
            nd=dict(dp)
            for s,c in dp.items(): nd[s+v]=nd.get(s+v,0)+c
            dp=nd
        w2=int(round(2*W)); dev=abs(w2-S/2)
        cnt=sum(c for s,c in dp.items() if abs(s-S/2)>=dev-1e-9)
        return cnt/N
    mu=tot/2; var=n*(n+1)*(2*n+1)/24-sum(t**3-t for t in ties)/48
    z=(abs(W-mu)-0.5)/math.sqrt(var) if var>0 else 0
    return math.erfc(max(z,0)/math.sqrt(2))
def test_vs_base(base,others,names):
    """base: list per unit; others: dict name->list. returns dict name->(adjp) and friedman p"""
    k=len(others)
    M=[[base[i]]+[others[nm][i] for nm in names] for i in range(len(base))]
    fp=friedman(M)
    res={}
    for nm in names:
        p=min(1.0,wilcoxon(others[nm],base)*k); res[nm]=p
    return fp,res
def mean_sd(v):
    m=sum(v)/len(v); s=math.sqrt(sum((a-m)**2 for a in v)/(len(v)-1)) if len(v)>1 else 0; return m,s
def load_dir(d,nch=None):
    out={}
    for f in sorted(glob.glob(RD+d+'/*.csv')):
        for r in csv.DictReader(open(f)):
            if nch is not None and r.get('Nch')!=str(nch): continue
            out.setdefault(r['Method'],{})[r['Rec']]=r
    return out
def load_one(fn):
    out={}
    for r in csv.DictReader(open(RD+fn)): out.setdefault(r['Method'],{})[r['Rec']]=r
    return out
