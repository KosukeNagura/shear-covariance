"""Independent mathematical and convention tests; no production simulations."""
import itertools
import numpy as np
import pytest
import shear_core as c


def fixture():
    ra=np.deg2rad([10,13,11,16]);dec=np.deg2rad([-8,-5,-2,1])
    v=np.column_stack((np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)))
    pairs=[np.array([[0,1],[0,2]]),np.array([[1,2],[2,3]])]
    w=np.array([1.,2.,.7,1.3])
    rng=np.random.default_rng(20261002);a=rng.normal(size=(8,8));s=a@a.T/8
    S=s.reshape(4,2,4,2).transpose(0,2,1,3)
    geo=[c.pair_geometry(p,v,w) for p in pairs]
    return v,pairs,w,s,S,geo


def test_matrix_blocks_against_quadratic_trace_and_direct_wick():
    v,pairs,w,s,S,geo=fixture()
    qs=[]
    for g in geo:
        per=[]
        for (i,j),weight,ri,rj in zip(*g):
            forms=[]
            for sign in [1,-1]:
                A=ri.T@np.diag([1,sign])@rj
                B=np.zeros((8,8));B[2*i:2*i+2,2*j:2*j+2]=weight*A/2
                B[2*j:2*j+2,2*i:2*i+2]=weight*A.T/2
                forms.append(B)
            per.append(forms)
        qs.append(np.sum(per,axis=0))
    for a in range(2):
        for b in range(2):
            expected=np.array([[2*np.trace(qs[a][x]@s@qs[b][y]@s) for y in range(2)] for x in range(2)])
            for lc,rc in [(0,2000),(1,1),(2,1)]:
                np.testing.assert_allclose(c.raw_pair_sum(geo[a],geo[b],S,rc,lc),expected,rtol=1e-12,atol=1e-12)
            if a==b:
                diagonal=np.einsum('i,i,iixy->xy',geo[a][1],geo[a][1],c.kernel_block(geo[a],geo[a],S))
                np.testing.assert_allclose(c.identical_sum(geo[a],S),diagonal,rtol=1e-12,atol=1e-12)
    # Reversing an unordered edge does not change its observable.
    rev=c.pair_geometry(pairs[0][:,::-1],v,w)
    np.testing.assert_allclose(c.raw_pair_sum(rev,geo[1],S),c.raw_pair_sum(geo[0],geo[1],S),atol=1e-12)


def test_ht_exact_expectation_all_subsets():
    v,pairs,w,s,S,geo=fixture();p=.001
    full=np.array([[c.raw_pair_sum(ga,gb,S) for gb in geo] for ga in geo])
    expectation=np.zeros_like(full)
    for bits in itertools.product([False,True],repeat=4):
        prob=p**sum(bits)*(1-p)**(4-sum(bits))
        sampled=[c.pair_geometry(pairs[b][np.array(bits[2*b:2*b+2])],v,w) for b in range(2)]
        for a in range(2):
            for b in range(2):
                raw=c.raw_pair_sum(sampled[a],sampled[b],S)
                identity=c.identical_sum(sampled[a],S) if a==b else 0
                expectation[a,b]+=prob*(raw-(1-p)*identity)/p**2
    np.testing.assert_allclose(expectation,full,rtol=3e-13,atol=1e-11)


def test_spherical_harmonic_addition_and_treecorr_convention():
    hp=pytest.importorskip('healpy');tc=pytest.importorskip('treecorr')
    ns=8;lm=4;pix=np.array([110,207,330,404]);v=np.array(hp.pix2vec(ns,pix)).T
    cl=np.array([0.,0.,.7,.4,.2]);S=np.zeros((8,8))
    # Sum a complete independent real E-harmonic basis: exact, no Monte Carlo tolerance.
    for ell in range(2,lm+1):
        for m in range(ell+1):
            for phase in ([1.] if m==0 else [1.,1j]):
                alm=np.zeros(hp.Alm.getsize(lm),complex)
                alm[hp.Alm.getidx(lm,ell,m)]=phase*np.sqrt(cl[ell] if m==0 else cl[ell]/2)
                maps=hp.alm2map([alm*0,alm,alm*0],ns,lmax=lm,pol=True)
                g=np.column_stack((maps[1][pix],-maps[2][pix])).ravel()
                S+=np.outer(g,g)
    th=np.linspace(0,np.pi,60001);xp,xm=c.xi_from_cl(th,cl)
    pred=c.signal_covariance(v,th,xp,xm,xp[0]).transpose(0,2,1,3).reshape(8,8)
    np.testing.assert_allclose(pred,S,atol=3e-10,rtol=3e-8)
    pairs=np.array(list(itertools.combinations(range(4),2)))
    g=np.array([[.3,.5],[-.7,.2],[.1,-.9],[.6,.8]]);w=np.array([1.,2.,3.,4.])
    geometry=c.pair_geometry(pairs,v,w)
    direct=c.measure_grid(g,[geometry],[geometry[1].sum()])
    ra,dec=hp.pix2ang(ns,pix,lonlat=True)
    gg=tc.GGCorrelation(min_sep=1,max_sep=179,nbins=1,sep_units='deg',metric='Arc',bin_slop=0,angle_slop=0)
    gg.process(tc.Catalog(ra=ra,dec=dec,ra_units='deg',dec_units='deg',w=w,g1=g[:,0],g2=g[:,1]),num_threads=1)
    np.testing.assert_allclose(direct,np.r_[gg.xip,gg.xim],atol=3e-8,rtol=3e-7)


def test_jacobi_recurrence_against_scipy_and_zero_lag():
    from scipy.special import eval_jacobi
    th=np.linspace(0,np.pi,71);cl=np.zeros(41);cl[2:]=1/(np.arange(2,41)+1)**2
    xp,xm=c.xi_from_cl(th,cl)
    for sign,answer,ab,factor in [(1,xp,(0,4),np.cos(th/2)**4),(-1,xm,(4,0),np.sin(th/2)**4)]:
        exact=sum((2*l+1)*cl[l]*eval_jacobi(l-2,*ab,np.cos(th)) for l in range(2,len(cl)))*factor/(4*np.pi)
        np.testing.assert_allclose(answer,exact,atol=1e-13,rtol=1e-11)
    assert xm[0]==0
    np.testing.assert_allclose(xp[0],sum((2*l+1)*cl[l] for l in range(2,len(cl)))/(4*np.pi))


def test_explicit_four_F_expressions():
    # Independent random D matrices and relative rotations; verify all signs.
    rng=np.random.default_rng(3)
    for _ in range(20):
        a,b,c0,d=rng.normal(size=4)
        p,m,q,n=rng.normal(size=4)
        def L(t):return np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
        h=L(a)@np.diag([(p+m)/2,(p-m)/2])@L(b).T
        h2=L(c0)@np.diag([(q+n)/2,(q-n)/2])@L(d).T
        de,se,df,sf=a-b,a+b,c0-d,c0+d
        F=.5*np.array([[p*q*np.cos(de-df)+m*n*np.cos(se-sf),p*n*np.cos(de-sf)+m*q*np.cos(se-df)],
                       [p*n*np.cos(de+sf)+m*q*np.cos(se+df),p*q*np.cos(de+df)+m*n*np.cos(se+sf)]])
        np.testing.assert_allclose(c.T@(h*h2)@c.T.T,F,atol=1e-14)
