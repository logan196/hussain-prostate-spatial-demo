/* Novaflow spatial stats engine (v2). Pure functions, no DOM. Runs in the browser and in Node.
   Every test here is computed from the served per-spot data; nothing is estimated by the language model.
   Conventions mirror squidpy 1.8 where a squidpy function exists:
   - nhood_enrichment: edge counts between label pairs on a radius graph, label-permutation null, z = (obs - mean) / sd (sd with ddof 0)
   - co_occurrence: squidpy's ratio p(A | B within d) / p(A within d), cumulative distance thresholds, self pairs excluded
   Permutations use a seeded mulberry32 stream with a Fisher-Yates shuffle, so a run is reproducible from (seed, n_perms). */

export function mulberry32(seed){let a=seed>>>0;return function(){a=(a+0x6D2B79F5)>>>0;let t=a;t=Math.imul(t^(t>>>15),t|1);t^=t+Math.imul(t^(t>>>7),t|61);return ((t^(t>>>14))>>>0)/4294967296}}
export function shuffle(arr,rnd){for(let i=arr.length-1;i>0;i--){const j=Math.floor(rnd()*(i+1));const t=arr[i];arr[i]=arr[j];arr[j]=t}return arr}

/* symmetric radius graph (CSR), no self loops. xs, ys in micrometres */
export function radiusGraph(xs,ys,r){const n=xs.length,cell=r,key=(cx,cy)=>cx*100003+cy,grid=new Map();
 for(let i=0;i<n;i++){const k=key(Math.floor(xs[i]/cell),Math.floor(ys[i]/cell));let b=grid.get(k);if(!b)grid.set(k,b=[]);b.push(i)}
 const r2=r*r,nbr=[];let m=0;
 for(let i=0;i<n;i++){const cx=Math.floor(xs[i]/cell),cy=Math.floor(ys[i]/cell),L=[];
  for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++){const b=grid.get(key(cx+dx,cy+dy));if(!b)continue;for(const j of b){if(j===i)continue;const ex=xs[i]-xs[j],ey=ys[i]-ys[j];if(ex*ex+ey*ey<=r2)L.push(j)}}
  L.sort((a,b)=>a-b);nbr.push(L);m+=L.length}
 const indptr=new Int32Array(n+1),indices=new Int32Array(m);let o=0;for(let i=0;i<n;i++){indptr[i]=o;for(const j of nbr[i])indices[o++]=j}indptr[n]=o;
 return {n,indptr,indices,edges:m,radius:r}}

function nhoodCount(g,lab,k){const c=new Float64Array(k*k);for(let i=0;i<g.n;i++){const a=lab[i]*k;for(let p=g.indptr[i];p<g.indptr[i+1];p++)c[a+lab[g.indices[p]]]++}return c}

/* labels: Int32Array of category codes 0..k-1 for every spot. Returns full k x k z matrix plus the null for pair (a,b). */
export function nhoodEnrichment(g,labels,k,{nPerms=1000,seed=0,pair=[0,1]}={}){
 const obs=nhoodCount(g,labels,k),sum=new Float64Array(k*k),sq=new Float64Array(k*k),pi=pair[0]*k+pair[1],nullAB=new Float64Array(nPerms);
 const rnd=mulberry32(seed),lab=Int32Array.from(labels);
 for(let t=0;t<nPerms;t++){shuffle(lab,rnd);const c=nhoodCount(g,lab,k);for(let q=0;q<k*k;q++){sum[q]+=c[q];sq[q]+=c[q]*c[q]}nullAB[t]=c[pi]}
 const z=new Float64Array(k*k),mean=new Float64Array(k*k),sd=new Float64Array(k*k);
 for(let q=0;q<k*k;q++){mean[q]=sum[q]/nPerms;sd[q]=Math.sqrt(Math.max(0,sq[q]/nPerms-mean[q]*mean[q]));z[q]=(obs[q]-mean[q])/sd[q]}
 let ge=0,le=0;for(const v of nullAB){if(v>=obs[pi])ge++;if(v<=obs[pi])le++}
 return {k,count:Array.from(obs),z:Array.from(z),mean:Array.from(mean),sd:Array.from(sd),pair,obs:obs[pi],zPair:z[pi],nullMean:mean[pi],nullSd:sd[pi],
  pEnriched:(ge+1)/(nPerms+1),pDepleted:(le+1)/(nPerms+1),nPerms,seed,null:Array.from(nullAB)}}

function median(a){const s=Float64Array.from(a).sort();const n=s.length;return n?(n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2):NaN}
/* nearest distance from every A spot to the closest B spot (self excluded), simple and exact */
export function nearestDists(xs,ys,aIdx,bIdx){const out=new Float64Array(aIdx.length),bx=Float64Array.from(bIdx,i=>xs[i]),by=Float64Array.from(bIdx,i=>ys[i]);
 for(let p=0;p<aIdx.length;p++){const i=aIdx[p],x=xs[i],y=ys[i];let best=Infinity;for(let q=0;q<bIdx.length;q++){if(bIdx[q]===i)continue;const dx=x-bx[q],dy=y-by[q],d=dx*dx+dy*dy;if(d<best)best=d}out[p]=Math.sqrt(best)}return out}

/* Proximity: is A closer to B than expected? statistic = mean nearest A->B distance (um); the median is reported too.
   Null: shuffle the label vector (scope "all": every spot in the section, squidpy style; scope "ab": only among A and B spots,
   which keeps the tissue footprint of the two groups fixed). One-sided p for "closer" and for "farther". */
export function proximity(xs,ys,labels,codeA,codeB,{nPerms=999,seed=0,withinUm=250,scope="all"}={}){
 const idx=(lab,c)=>{const o=[];for(let i=0;i<lab.length;i++)if(lab[i]===c)o.push(i);return o};
 const mean=a=>{let s=0;for(const v of a)s+=v;return s/a.length};
 const A=idx(labels,codeA),B=idx(labels,codeB),d=nearestDists(xs,ys,A,B),obs=mean(d);let within=0;for(const v of d)if(v<=withinUm)within++;
 const rnd=mulberry32(seed),lab=Int32Array.from(labels),nul=new Float64Array(nPerms);
 const pos=[];if(scope==="ab")for(let i=0;i<lab.length;i++)if(lab[i]===codeA||lab[i]===codeB)pos.push(i);const sub=Int32Array.from(pos,i=>labels[i]);
 for(let t=0;t<nPerms;t++){if(scope==="ab"){shuffle(sub,rnd);for(let q=0;q<pos.length;q++)lab[pos[q]]=sub[q]}else shuffle(lab,rnd);nul[t]=mean(nearestDists(xs,ys,idx(lab,codeA),idx(lab,codeB)))}
 let le=0,ge=0,s=0,s2=0;for(const v of nul){if(v<=obs)le++;if(v>=obs)ge++;s+=v;s2+=v*v}const m=s/nPerms,sd=Math.sqrt(Math.max(0,s2/nPerms-m*m));
 return {nA:A.length,nB:B.length,obs,median_obs:median(d),fracWithin:within/d.length,withinUm,scope,nullMean:m,nullSd:sd,z:(obs-m)/sd,
  pCloser:(le+1)/(nPerms+1),pFarther:(ge+1)/(nPerms+1),nPerms,seed,null:Array.from(nul),dists:Array.from(d)}}

/* squidpy co_occurrence, exact: counts over ordered pairs i!=j with d^2 <= t^2 (cumulative), occ[a,b] = (cnt[b,a]/rowsum[b]) / (rowsum[a]/total) */
export function coOccurrence(xs,ys,labels,k,thresholds){const L=thresholds.length,n=xs.length,cnt=new Float64Array(k*k*L),t2=thresholds.map(t=>t*t);
 for(let i=0;i<n;i++){const xi=xs[i],yi=ys[i],li=labels[i]*k;for(let j=0;j<n;j++){if(i===j)continue;const dx=xi-xs[j],dy=yi-ys[j],d2=dx*dx+dy*dy;if(d2>t2[L-1])continue;const base=(li+labels[j])*L;for(let r=0;r<L;r++)if(d2<=t2[r])cnt[base+r]++}}
 return occFromCounts(cnt,k,L)}
function occFromCounts(cnt,k,L){const occ=new Float64Array(k*k*L);
 for(let r=0;r<L;r++){const row=new Float64Array(k);let tot=0;for(let c=0;c<k;c++){let s=0;for(let i=0;i<k;i++)s+=cnt[(i*k+c)*L+r];row[c]=s;tot+=s}
  for(let c=0;c<k;c++)for(let i=0;i<k;i++){const pr=row[i]/tot;if(pr!==0&&row[c]!==0)occ[(i*k+c)*L+r]=(cnt[(c*k+i)*L+r]/row[c])/pr}}
 return {k,L,occ:Array.from(occ),at:(a,b,r)=>occ[(a*k+b)*L+r]}}

/* co-occurrence of A given B at one distance band, with a label-permutation null over a precomputed pair list */
export function coOccurrenceTest(xs,ys,labels,k,codeA,codeB,thresholds,{distanceUm,nPerms=499,seed=0}={}){
 const full=coOccurrence(xs,ys,labels,k,thresholds),L=thresholds.length,curve=thresholds.map((t,r)=>full.at(codeA,codeB,r));
 const r0=distanceUm,r2=r0*r0,n=xs.length,I=[],J=[];
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){if(i===j)continue;const dx=xs[i]-xs[j],dy=ys[i]-ys[j];if(dx*dx+dy*dy<=r2){I.push(i);J.push(j)}}
 const PI=Int32Array.from(I),PJ=Int32Array.from(J);
 const stat=lab=>{const c=new Float64Array(k*k);for(let p=0;p<PI.length;p++)c[lab[PI[p]]*k+lab[PJ[p]]]++;const one=occFromCounts(Float64Array.from(c),k,1);return one.at(codeA,codeB,0)};
 const obs=stat(labels),rnd=mulberry32(seed),lab=Int32Array.from(labels),nul=new Float64Array(nPerms);
 for(let t=0;t<nPerms;t++){shuffle(lab,rnd);nul[t]=stat(lab)}
 let ge=0,le=0,s=0,s2=0;for(const v of nul){if(v>=obs)ge++;if(v<=obs)le++;s+=v;s2+=v*v}const mean=s/nPerms,sd=Math.sqrt(Math.max(0,s2/nPerms-mean*mean));
 return {thresholds,curve,distanceUm:r0,obs,nullMean:mean,nullSd:sd,z:(obs-mean)/sd,pHigher:(ge+1)/(nPerms+1),pLower:(le+1)/(nPerms+1),nPerms,seed,pairs:PI.length,null:Array.from(nul)}}

/* normal upper tail, accurate in the far tail; returns {p, log10p} */
export function normSf(z){z=Math.abs(z);const x=z/Math.SQRT2;
 if(z<37){const p=0.5*erfc(x);return {p,log10p:Math.log10(p)}}
 let f=0;for(let n=80;n>=1;n--)f=(n/2)/(x+f);const lnp=Math.log(0.5)-x*x-0.5*Math.log(Math.PI)-Math.log(x+f);return {p:Math.exp(lnp),log10p:lnp/Math.LN10}}
function erfc(x){/* x >= 0. Series for erf below 2, continued fraction above; relative error well under 1e-10 */
 if(x<2){let term=x,sum=x;for(let n=1;n<80;n++){term*=-x*x/n;const add=term/(2*n+1);sum+=add;if(Math.abs(add)<1e-17*Math.abs(sum))break}return 1-2/Math.sqrt(Math.PI)*sum}
 let f=0;for(let n=80;n>=1;n--)f=(n/2)/(x+f);return Math.exp(-x*x)/Math.sqrt(Math.PI)/(x+f)}

/* Wilcoxon rank-sum (Mann-Whitney U), two-sided, normal approximation with tie correction and continuity correction (scipy.stats.mannwhitneyu defaults, method="asymptotic") */
export function wilcoxon(a,b){const n1=a.length,n2=b.length,N=n1+n2,all=new Float64Array(N);all.set(a,0);all.set(b,n1);
 const ord=Array.from({length:N},(_,i)=>i).sort((i,j)=>all[i]-all[j]),rank=new Float64Array(N);let tie=0;
 for(let i=0;i<N;){let j=i;while(j+1<N&&all[ord[j+1]]===all[ord[i]])j++;const r=(i+j+2)/2;for(let q=i;q<=j;q++)rank[ord[q]]=r;const t=j-i+1;tie+=t*t*t-t;i=j+1}
 let R1=0;for(let i=0;i<n1;i++)R1+=rank[i];const U1=R1-n1*(n1+1)/2,mu=n1*n2/2,s=Math.sqrt(n1*n2/12*((N+1)-tie/(N*(N-1))));
 const U=Math.max(U1,n1*n2-U1),z=(U-mu-0.5)/s,sf=normSf(z),p=Math.min(1,2*sf.p),log10p=Math.min(0,Math.log10(2)+sf.log10p);
 const mean=v=>v.reduce((x,y)=>x+y,0)/v.length,pos=v=>v.filter(x=>x>0).length/v.length,em=v=>v.reduce((x,y)=>x+Math.expm1(y),0)/v.length;
 const q=(v,f)=>{const s=Float64Array.from(v).sort();const h=(s.length-1)*f,l=Math.floor(h);return s[l]+(s[Math.min(l+1,s.length-1)]-s[l])*(h-l)};
 return {n1,n2,U1,z:(U1-mu)/s,zc:z,p,log10p,auc:U1/(n1*n2),meanA:mean(a),meanB:mean(b),pctA:pos(a),pctB:pos(b),log2fc:Math.log2((em(a)+1e-9)/(em(b)+1e-9)),
  qA:[0.1,0.25,0.5,0.75,0.9].map(f=>q(a,f)),qB:[0.1,0.25,0.5,0.75,0.9].map(f=>q(b,f))}}
