/* v2: live LLM chat that drives the Vitessce view + hypothesis tests computed in code.
   The model (server-side, /api/ask) returns a structured "respond" tool call: view changes and test choices.
   Every number on a result card comes from v2/stats.js running in a Web Worker over the served data. */
const ENDPOINT=(document.querySelector('meta[name="llm-endpoint"]')||{}).content||"api/ask";
const COL={a:[236,72,153],b:[79,140,255],c:[245,166,35],d:[52,211,153]};
const hex=c=>"#"+c.map(v=>v.toString(16).padStart(2,"0")).join("");
const rgb=h=>{const m=/^#?([0-9a-f]{6})$/i.exec(h||"");return m?[0,2,4].map(i=>parseInt(m[1].slice(i,i+2),16)):null};
const esc=s=>String(s==null?"":s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const f2=(v,d=2)=>v==null||!isFinite(v)?"n/a":(+v).toFixed(d);
function fp(p,log10p){if(p==null)return "n/a";if(p>=0.001)return p.toPrecision(2).replace(/\.?0+$/,"");if(p>0)return p.toExponential(1);return "< 1e"+Math.floor(log10p||-300)}
let A=null,N=0,XS=null,YS=null,worker=null,wid=0;const pending=new Map();
const turns=[];let health=null;window.__results=[];

/* ---------------- data ---------------- */
let bloscP=null;const chunks=new Map();
async function geneValues(sym){await A.genesReady;const gi=A.sym(sym);if(gi==null)return null;const ch=Math.floor(gi/40),off=gi%40;
 if(!chunks.has(ch))chunks.set(ch,(async()=>{if(!bloscP)bloscP=import("./vendor/blosc.js").then(m=>m.default.fromConfig({id:"blosc",cname:"zstd",clevel:9,shuffle:2,blocksize:0}));
  const codec=await bloscP,r=await fetch(A.base+"vitessce/data.zarr/X/0."+ch);if(!r.ok)throw new Error("matrix chunk "+ch+" HTTP "+r.status);
  const out=await codec.decode(new Uint8Array(await r.arrayBuffer()));return new Float32Array(out.buffer,out.byteOffset,out.byteLength/4)})());
 const f=await chunks.get(ch),v=new Float64Array(N);for(let i=0;i<N;i++)v[i]=f[i*40+off];return v}
async function resolveGene(s){if(!s)return null;await A.genesReady;const u=String(s).trim().toUpperCase();if(A.sym(u)!=null)return {g:A.rowSym(A.sym(u)),how:"exact"};
 const al=A.alias(u);if(al)return {g:al,how:"alias "+s};const f=A.findGene(String(s));return f||null}
const STATES=["AR adenocarcinoma","NE / t-SCNC-like","Amphicrine (AR+NE)","Stroma","Other"];
const PROG={"ar-high":0,"ne-high":5,"dual ar+/ne+":2};
function normState(v){const t=String(v||"").toLowerCase();if(/abstain|uncertain/.test(t))return -1;if(/amphi|dual|hybrid/.test(t))return 2;if(/strom/.test(t))return 3;if(/^ne|neuro|scnc|small/.test(t))return 1;if(/^ar|adeno|androgen/.test(t))return 0;if(/other/.test(t))return 4;return null}
async function resolveGroup(spec,other){if(!spec||!spec.kind)throw new Error("missing group");const m=new Uint8Array(N);let label=spec.label||"",how="";const k=spec.kind,op=spec.op||">",th=+(spec.threshold??0);
 const cmp=(v)=>op===">"?v>th:op===">="?v>=th:op==="<"?v<th:v<=th;
 if(k==="state"){const c=normState(spec.value);if(c==null)throw new Error(`unknown state "${spec.value}"`);for(let i=0;i<N;i++)m[i]=A.SP.state[i]===c;label=label||(c<0?"Abstained":STATES[c]);how="ensemble state"}
 else if(k==="program"){const t=String(spec.value||"").toLowerCase().replace(/\s+/g,"");let code=null;if(/^ar/.test(t))code=0;else if(/^ne/.test(t))code=5;else if(/dual/.test(t))code=2;
  if(code!=null)for(let i=0;i<N;i++)m[i]=A.SP.region[i]===code;else if(/low|other/.test(t))for(let i=0;i<N;i++)m[i]=![0,5,2].includes(A.SP.region[i]);else throw new Error(`unknown program "${spec.value}"`);
  label=label||({0:"AR-high",5:"NE-high",2:"Dual AR+/NE+"}[code]||"Low/other");how="program z>1"}
 else if(k==="region"){const idx=A.MT.regions.findIndex(r=>r.toLowerCase()===String(spec.value||"").toLowerCase());if(idx<0)throw new Error(`unknown region "${spec.value}"`);for(let i=0;i<N;i++)m[i]=A.SP.region[i]===idx;label=label||A.MT.regions[idx];how="program region"}
 else if(k==="leiden"){const c=String(spec.value).replace(/\D/g,"");for(let i=0;i<N;i++)m[i]=String(A.MT.leiden[i])===c;label=label||"Leiden "+c;how="Leiden cluster"}
 else if(k==="gene"){const g=await resolveGene(spec.gene||spec.value);if(!g)throw new Error(`gene "${spec.gene||spec.value}" is not in the 18,939-gene matrix`);const v=await geneValues(g.g);for(let i=0;i<N;i++)m[i]=cmp(v[i]);label=label||`${g.g} ${op} ${th}`;how=`${g.g} log-norm expression ${op} ${th}`+(g.how!="exact"?` (${g.how})`:"")}
 else if(k==="score"){const v=String(spec.score||spec.value).toUpperCase().startsWith("N")?A.MT.NE_z:A.MT.AR_z;const nm=v===A.MT.NE_z?"NE":"AR";for(let i=0;i<N;i++)m[i]=cmp(v[i]);label=label||`${nm} score ${op} ${th}`;how=`${nm} program z ${op} ${th}`}
 else if(k==="distance_to_ar"){const o=spec.op||"<=";if(!DAR){const ar=[];for(let i=0;i<N;i++)if(A.SP.region[i]===0)ar.push(i);DAR=new Float64Array(N);for(let i=0;i<N;i++){let b=Infinity;for(const j of ar){const dx=XS[i]-XS[j],dy=YS[i]-YS[j],q=dx*dx+dy*dy;if(q<b)b=q}DAR[i]=Math.sqrt(b)}}for(let i=0;i<N;i++){const d=DAR[i];m[i]=o==="<="?d<=th:o==="<"?d<th:o===">"?d>th:d>=th}label=label||`${o} ${th} µm from AR`;how="distance to nearest AR-high spot"}
 else if(k==="rest"){if(!other)throw new Error("'rest' needs another group");for(let i=0;i<N;i++)m[i]=!other[i];label=label||"all other spots";how="complement"}
 else throw new Error(`unknown group kind "${k}"`);
 if(Array.isArray(spec.and))for(const s of spec.and){const r=await resolveGroup(s,other);for(let i=0;i<N;i++)m[i]=m[i]&&r.mask[i];label=spec.label||label+" ∧ "+r.label;how+=" ∧ "+r.how}
 let n=0;for(let i=0;i<N;i++)n+=m[i];return {mask:m,label,how,n}}
function labelsFor(ma,mb){const lab=new Int32Array(N);let rest=0,both=0;for(let i=0;i<N;i++){if(ma[i]&&mb[i]){lab[i]=-2;both++}else if(ma[i])lab[i]=0;else if(mb[i])lab[i]=1;else{lab[i]=-1;rest++}}
 const cats=["A","B"];let rc=-1,bc=-1;if(rest){rc=cats.length;cats.push("rest")}if(both){bc=cats.length;cats.push("A and B")}
 for(let i=0;i<N;i++){if(lab[i]===-1)lab[i]=rc;else if(lab[i]===-2)lab[i]=bc}return {labels:lab,k:cats.length,cats,both,rest}}
function run(msg,transfer){return new Promise((res,rej)=>{const id=++wid;pending.set(id,{res,rej});worker.postMessage({...msg,id},transfer||[])})}

/* ---------------- tests ---------------- */
const DEF={nhood_enrichment:{n:1000,radius:130},proximity:{n:999,within:250},co_occurrence:{n:499,dist:300},wilcoxon:{}};
const THR=[150,300,450,600,800,1000,1500,2000];
async function runTest(t){const kind=t.test,id="t"+(window.__results.length+1);
 const ga=await resolveGroup(t.group_a,null),gb=await resolveGroup(t.group_b,ga.mask);
 if(!ga.n||!gb.n)throw new Error(`empty group: ${ga.label} n=${ga.n}, ${gb.label} n=${gb.n}`);
 const {labels,k,cats,both}=labelsFor(ga.mask,gb.mask);const seed=0;let out,params={},gene=null;
 if(kind==="wilcoxon"){const g=await resolveGene(t.gene);if(!g)throw new Error(`gene "${t.gene}" is not in the matrix`);gene=g.g;const v=await geneValues(gene),a=[],b=[];
  for(let i=0;i<N;i++){if(labels[i]===0)a.push(v[i]);else if(labels[i]===1)b.push(v[i])}out=await run({test:kind,a,b});params={gene,gene_how:g.how}}
 else if(kind==="nhood_enrichment"){const r=Math.min(400,Math.max(60,+t.radius_um||DEF[kind].radius)),n=Math.min(2000,Math.max(99,+t.n_perms||DEF[kind].n));params={radius_um:r,n_perms:n,seed};
  out=await run({test:kind,xs:XS,ys:YS,labels,k,radius:r,nPerms:n,seed})}
 else if(kind==="proximity"){const n=Math.min(2000,Math.max(99,+t.n_perms||DEF[kind].n)),scope=t.null_scope==="ab"?"ab":"all";params={n_perms:n,seed,within_um:DEF[kind].within,null_scope:scope};out=await run({test:kind,xs:XS,ys:YS,labels,k,nPerms:n,seed,withinUm:DEF[kind].within,scope})}
 else if(kind==="co_occurrence"){const n=Math.min(2000,Math.max(99,+t.n_perms||DEF[kind].n)),d=Math.min(1000,Math.max(120,+t.distance_um||DEF[kind].dist));params={n_perms:n,seed,distance_um:d,interval_um:[0,...THR]};
  out=await run({test:kind,xs:XS,ys:YS,labels,k,thresholds:THR,distanceUm:d,nPerms:n,seed})}
 else throw new Error("unknown test "+kind);
 const idxA=[],idxB=[];for(let i=0;i<N;i++){if(ga.mask[i])idxA.push(i);if(gb.mask[i])idxB.push(i)}
 let uA=0,uB=0;for(let i=0;i<N;i++){if(labels[i]===0)uA++;else if(labels[i]===1)uB++}
 const rec={id,test:kind,a:{label:ga.label,how:ga.how,n:ga.n,n_used:uA},b:{label:gb.label,how:gb.how,n:gb.n,n_used:uB},overlap:both,cats,params,gene,result:out,idxA,idxB,spec:t};
 window.__results.push(rec);return rec}

/* ---------------- result cards ---------------- */
function hist(nul,obs,{w=300,h=86,xl=""}={}){if(!nul||!nul.length)return "";let lo=Math.min(...nul),hi=Math.max(...nul);if(hi===lo){hi+=1;lo-=1}
 if(lo>=0&&obs>=0&&nul.every(v=>v>=0))nul.nonneg=true;const span=hi-lo,far=obs>hi+4*span?1:obs<lo-4*span?-1:0;if(!far){lo=Math.min(lo,obs);hi=Math.max(hi,obs)}const pad=(hi-lo)*.06;lo-=pad;hi+=pad;if(nul.nonneg)lo=Math.max(0,lo);
 const L=far<0?62:6,R=far>0?62:6,T=6,Bt=22,B=28,c=new Array(B).fill(0);for(const v of nul)c[Math.min(B-1,Math.max(0,Math.floor((v-lo)/(hi-lo)*B)))]++;
 const mx=Math.max(...c),X=v=>L+(v-lo)/(hi-lo)*(w-L-R),bw=(w-L-R)/B,d=v=>f2(v,Number.isInteger(v)||Math.abs(v)>=10?0:2);
 const ox=far>0?w-8:far<0?8:X(obs),brk=far?`<text x="${far>0?w-R+14:L-14}" y="${(h-Bt)/2+4}" text-anchor="middle" fill="#9aa6b6">⋯</text>`:"";
 return `<svg viewBox="0 0 ${w} ${h}" width="100%" class="mini">${c.map((n,i)=>`<rect x="${(L+i*bw).toFixed(1)}" y="${(T+(1-n/mx)*(h-T-Bt)).toFixed(1)}" width="${Math.max(.5,bw-1).toFixed(1)}" height="${(n/mx*(h-T-Bt)).toFixed(1)}" fill="#4a5668"/>`).join("")}${brk}
 <line x1="${ox}" x2="${ox}" y1="${T-2}" y2="${h-Bt}" stroke="#ec4899" stroke-width="2"/><text x="${ox+(far>0?-4:far<0?4:0)}" y="${far?T+8:h-Bt+11}" text-anchor="${far>0||(!far&&ox>w*.7)?"end":far<0||ox<w*.3?"start":"middle"}" fill="#ec4899">observed ${d(obs)}</text>
 <text x="${L}" y="${h-2}">${d(lo)}</text><text x="${w-R}" y="${h-2}" text-anchor="end">${d(hi)}</text><text x="${(L+w-R)/2}" y="${h-2}" text-anchor="middle">${esc(xl)}</text></svg>`}
function curve(th,ys,mark,{w=300,h=92}={}){const L=26,R=6,T=6,Bt=22,xmax=th[th.length-1],vals=ys.filter(isFinite),hi=Math.max(1.2,...vals)*1.08,lo=0,X=v=>L+v/xmax*(w-L-R),Y=v=>T+(1-(v-lo)/(hi-lo))*(h-T-Bt);
 const p=th.map((t,i)=>`${i?"L":"M"}${X(t).toFixed(1)},${Y(isFinite(ys[i])?ys[i]:0).toFixed(1)}`).join("");
 return `<svg viewBox="0 0 ${w} ${h}" width="100%" class="mini"><line x1="${L}" x2="${w-R}" y1="${Y(1)}" y2="${Y(1)}" stroke="#555" stroke-dasharray="3 3"/><text x="${L-4}" y="${Y(1)+3}" text-anchor="end">1</text>
 <path d="${p}" stroke="#ec4899" stroke-width="2" fill="none"/>${th.map((t,i)=>`<circle cx="${X(t)}" cy="${Y(isFinite(ys[i])?ys[i]:0)}" r="2.4" fill="#ec4899"/>`).join("")}
 <line x1="${X(mark)}" x2="${X(mark)}" y1="${T}" y2="${h-Bt}" stroke="#38bdf8" stroke-dasharray="2 3"/><text x="${L}" y="${h-2}">0</text><text x="${w-R}" y="${h-2}" text-anchor="end">${xmax} µm</text><text x="${(L+w)/2}" y="${h-2}" text-anchor="middle">distance</text></svg>`}
function boxes(r,la,lb,{w=300,h=74}={}){const all=[...r.qA,...r.qB],hi=Math.max(...all,0.01),L=6,R=6,X=v=>L+v/hi*(w-L-R);
 const row=(q,y,c,l)=>`<line x1="${X(q[0])}" x2="${X(q[4])}" y1="${y}" y2="${y}" stroke="${c}"/><rect x="${X(q[1])}" y="${y-7}" width="${Math.max(1,X(q[3])-X(q[1]))}" height="14" fill="${c}" fill-opacity=".35" stroke="${c}"/><line x1="${X(q[2])}" x2="${X(q[2])}" y1="${y-8}" y2="${y+8}" stroke="#fff" stroke-width="2"/><text x="${L}" y="${y-10}">${esc(l)}</text>`;
 return `<svg viewBox="0 0 ${w} ${h}" width="100%" class="mini">${row(r.qA,20,"#ec4899",la)}${row(r.qB,50,"#4f8cff",lb)}<text x="${L}" y="${h-1}">0</text><text x="${w-R}" y="${h-1}" text-anchor="end">${f2(hi)} (log-norm)</text></svg>`}
const CIRC=/^(AR|KLK3|KLK2|TMPRSS2|NKX3-1|FKBP5|STEAP2|PMEPA1|SLC45A3|HOXB13|FOLH1|CHGA|CHGB|SYP|ASCL1|INSM1|NCAM1|ENO2|SOX2|SCG2|SCG3|PCSK1|STMN2|ELAVL4|NEUROD1)$/;
function card(rec){const r=rec.result,a=rec.a,b=rec.b,P=rec.params,rows=[];let title="",sent="",plot="",method="",cav=[];
 const n=`n<sub>A</sub> = ${a.n_used.toLocaleString()} · n<sub>B</sub> = ${b.n_used.toLocaleString()} spots`+(rec.overlap?` (${rec.overlap} spots in both groups set aside as their own category; A had ${a.n}, B had ${b.n})`:"");
 if(rec.test==="proximity"){title="Spatial proximity";const dir=r.pCloser<0.05?"closer to":r.pFarther<0.05?"farther from":null;
  sent=dir?`${a.label} spots sit ${dir} ${b.label} spots than label shuffling predicts.`:`No departure from the shuffled-label null for ${a.label} to ${b.label} distance.`;
  rows.push(["Statistic",`mean nearest A→B distance <b>${f2(r.obs,0)} µm</b> (median ${f2(r.median_obs,0)}; ${(100*r.fracWithin).toFixed(1)}% within ${r.withinUm} µm)`],["Null",`${r.scope==="ab"?"labels shuffled among A and B spots only":"labels shuffled across all spots"}: ${f2(r.nullMean,0)} ± ${f2(r.nullSd,1)} µm · z = <b>${f2(r.z,1)}</b>`],
   ["p",`closer: <b>${fp(r.pCloser)}</b> · farther: ${fp(r.pFarther)} (one-sided, (k+1)/(N+1))`],["n",n],["Permutations",`${r.nPerms} · seed ${r.seed} · ${r.ms} ms`]);
  plot=hist(r.null,r.obs,{xl:"null mean distance (µm)"});method="Mean distance from each A spot to its nearest B spot, compared with the same statistic after shuffling spot labels (group sizes kept). The smallest reachable p is 1/(N+1).";}
 else if(rec.test==="nhood_enrichment"){title="Neighbourhood enrichment";const dir=r.pEnriched<0.05?"more":r.pDepleted<0.05?"fewer":null;
  sent=dir?`${a.label} spots touch ${b.label} spots ${dir} often than label shuffling predicts.`:`A–B contacts are within the shuffled-label range.`;
  rows.push(["Statistic",`A–B neighbour pairs <b>${r.obs.toLocaleString()}</b> · z = <b>${f2(r.zPair)}</b>`],["Null",`label permutation: ${f2(r.nullMean,1)} ± ${f2(r.nullSd,1)} pairs`],["p",`enriched: <b>${fp(r.pEnriched)}</b> · depleted: ${fp(r.pDepleted)} (one-sided)`],
   ["n",n],["Graph",`radius ${P.radius_um} µm · ${r.edges.toLocaleString()} directed edges · mean degree ${f2(r.meanDegree,1)}`],["Permutations",`${r.nPerms} · seed ${r.seed} · ${r.ms} ms`]);
  plot=hist(r.null,r.obs,{xl:"null A–B neighbour pairs"});method="squidpy.gr.nhood_enrichment logic: count A–B edges on a radius graph, z against label-shuffled counts (sd with ddof 0).";}
 else if(rec.test==="co_occurrence"){title="Co-occurrence";const dir=r.pHigher<0.05?"above":r.pLower<0.05?"below":null;
  sent=dir?`Within ${r.distanceUm} µm of ${b.label}, ${a.label} is ${dir} its section-wide share.`:`Within ${r.distanceUm} µm of ${b.label}, ${a.label} is within the shuffled-label range.`;
  rows.push(["Statistic",`ratio p(A | within ${r.distanceUm} µm of B) / p(A) = <b>${f2(r.obs,3)}</b>`],["Null",`label permutation: ${f2(r.nullMean,3)} ± ${f2(r.nullSd,3)} · z = <b>${f2(r.z)}</b>`],["p",`higher: <b>${fp(r.pHigher)}</b> · lower: ${fp(r.pLower)} (one-sided)`],["n",n+` · ${r.pairs.toLocaleString()} ordered pairs ≤ ${r.distanceUm} µm`],["Permutations",`${r.nPerms} · seed ${r.seed} · ${r.ms} ms`]);
  plot=curve(r.thresholds,r.curve,r.distanceUm);method="squidpy.gr.co_occurrence ratio (cumulative distance thresholds, self pairs excluded); ratio 1 = same as section-wide share. Curve shows the ratio at each threshold.";}
 else if(rec.test==="wilcoxon"){title=`Differential expression · ${rec.gene}`;const dir=r.p<0.05?(r.auc>0.5?"higher":"lower"):null;
  sent=dir?`${rec.gene} is ${dir} in ${a.label} than in ${b.label} at the spot level.`:`No spot-level difference in ${rec.gene} between the two groups at α = 0.05.`;
  rows.push(["Statistic",`U = ${r.U1.toLocaleString()} · AUC = <b>${f2(r.auc,3)}</b> · log2FC = <b>${f2(r.log2fc)}</b>`],["Means",`A ${f2(r.meanA,3)} (${(100*r.pctA).toFixed(0)}% > 0) · B ${f2(r.meanB,3)} (${(100*r.pctB).toFixed(0)}% > 0)`],
   ["Null","no difference in distribution (rank-sum, normal approximation, tie + continuity corrected)"],["p",`two-sided <b>${fp(r.p,r.log10p)}</b> · z = ${f2(r.z)}`],["n",n]);
  plot=boxes(r,a.label,b.label);method="Wilcoxon rank-sum / Mann-Whitney U as scipy.stats.mannwhitneyu (asymptotic, two-sided).";
  cav.push("Spots are not independent replicates: neighbouring spots share tissue and signal, so this p-value overstates certainty. Treat it as a description of this section.");
  if(CIRC.test(rec.gene)&&/program|score|state|region/.test([rec.spec.group_a.kind,rec.spec.group_b.kind].join(" ")))cav.push(`${rec.gene} is one of the genes used to define the AR/NE programs, so comparing it across program-defined groups is partly circular.`)}
 cav.push("Unit is the Visium spot (55 µm, several cells), one section from one patient: no claim beyond this tissue, nothing causal or clinical.");
 if(rec.test!=="wilcoxon")cav.push("Distances use the array's own 100 µm spot pitch as the scale.");
 if(rec.test!=="wilcoxon")cav.push("The label-shuffle null treats every spot position as exchangeable; it does not preserve tissue architecture.");
 return `<div class="rc" data-id="${rec.id}"><div class="rch"><span>${title}</span><span class="tag live">COMPUTED</span></div>
 <div class="rcg"><i style="background:#ec4899"></i>A: ${esc(a.label)} <small>(${esc(a.how)})</small><br><i style="background:#4f8cff"></i>B: ${esc(b.label)} <small>(${esc(b.how)})</small></div>
 <p class="rcs">${esc(sent)}</p>${plot}<table class="rct">${rows.map(([k,v])=>`<tr><td>${k}</td><td>${v}</td></tr>`).join("")}</table>
 <div class="rcm"><b>Method.</b> ${esc(method)}</div><ul class="rcc">${cav.map(c=>`<li>${esc(c)}</li>`).join("")}</ul>
 <button class="rcx" data-id="${rec.id}">copy result JSON</button></div>`}
function pendingCard(t){return `<div class="rc"><div class="rch"><span>${esc(t.test.replace("_"," "))}</span><span class="tag">COMPUTING</span></div><div class="spin">running ${t.n_perms||DEF[t.test]?.n||""} ${t.test==="wilcoxon"?"":"permutations "}in your browser…</div></div>`}

/* ---------------- view actions ---------------- */
async function applyView(v,tests){const S=A.state,notes=[];v=v||{};
 if(v.preset&&A.LAYOUTS[v.preset]){A.preset(v.preset,false,true);notes.push("preset: "+v.preset)}
 let colored=false;const cb=v.color_by||{};
 if(cb.kind==="gene"&&cb.gene){const g=await resolveGene(cb.gene);if(g){S.enc="geneSelection";S.gene=g.g;S.extra=null;S.extraColors=null;colored=true;notes.push("colour: "+g.g+(g.how!="exact"?` (${g.how})`:""))}else notes.push(`unknown gene "${cb.gene}"`)}
 else if(cb.kind==="score"){const m={AR:"◆ AR score (z)",NE:"◆ NE score (z)",stroma:"◆ stroma score (z)",confidence:"◆ state confidence",distance_to_ar:"◆ distance to AR (mm)"}[cb.score||"NE"];S.enc="geneSelection";S.gene=m;S.extra=null;S.extraColors=null;colored=true;notes.push("colour: "+m.slice(2))}
 else if(cb.kind==="state"){S.enc="cellSetSelection";S.extra=null;S.extraColors=null;S.sel=Object.keys(A.STATECOL).map(k=>A.CS(k));colored=true;notes.push("colour: cell state")}
 else if(cb.kind==="region"){S.enc="cellSetSelection";S.extra=null;S.extraColors=null;S.sel=A.MT.regions.map(r=>[A.SETS[1],r]);colored=true;notes.push("colour: program region")}
 else if(cb.kind==="leiden"){S.enc="cellSetSelection";S.extra=null;S.extraColors=null;S.sel=A.MT.clusters.map(c=>[A.SETS[2],c]);colored=true;notes.push("colour: Leiden")}
 let hl=(v.highlight||[]).slice(0,4);const groups=[];
 if(tests&&tests.length){const t=tests[0];hl=[{group:t.group_a,color:"#ec4899"},{group:t.group_b,color:"#4f8cff"}]}
 let prev=null;for(const [i,h] of hl.entries()){try{const r=await resolveGroup(h.group,prev);prev=r.mask;groups.push({r,c:rgb(h.color)||Object.values(COL)[i%4]})}catch(e){notes.push(e.message)}}
 const geneWins=colored&&(cb.kind==="gene"||cb.kind==="score")&&(!tests||!tests.length||tests.some(t=>t.test==="wilcoxon"));
 if(groups.length&&!geneWins){if(colored)notes.push("highlight shown instead of the colour, to display the tested groups");const sets=groups.map(({r,c})=>[`${r.label} (${r.n})`,A.SP.barcode.filter((_,i)=>r.mask[i]),c]);A.setExtra("Ask",sets);S.sel=sets.map(([n])=>["Ask",n]);notes.push("highlight: "+groups.map(g=>`${g.r.label} (${g.r.n})`).join(", "))}
 const z=(v.zoom_to||"").toLowerCase();const cols=A.LAYOUTS[S.preset][0][3];
 if(z&&z!=="all"&&groups.length){const m=z==="b"&&groups[1]?groups[1].r.mask:z==="a"?groups[0].r.mask:null;const f=m?(i=>m[i]):(i=>groups.some(g=>g.r.mask[i]));S.view=A.fit(cols,f);notes.push("zoom: "+(m?(z==="a"?groups[0]:groups[1]).r.label:"highlighted spots"))}
 else if(z==="all"){S.view=A.fit(cols,null);notes.push("zoom: whole section")}
 A.render();return notes}

/* ---------------- chat panel ---------------- */
function panel(){const st=health==null?`<span class="llmpill">checking LLM…</span>`:health.key_configured?`<span class="llmpill on" title="${esc(health.model||"")}">live LLM</span>`:`<span class="llmpill off" title="ANTHROPIC_API_KEY is not set on this deployment">LLM key not set · rule parser</span>`;
 const T=turns.slice(-4).map(t=>`<div class="turn"><div class="q">“${esc(t.q)}”</div>${t.status==="thinking"?`<div class="spin">asking the model…</div>`:""}
  ${t.err?`<div class="err">${esc(t.err)}</div>`:""}${t.answer?`<p class="ans">${esc(t.answer)}</p>`:""}${t.interp?`<p class="interp">${esc(t.interp)}</p>`:""}
  ${t.notes&&t.notes.length?`<div class="pchips">${t.notes.map(n=>`<span>${esc(n)}</span>`).join("")}</div>`:""}${(t.cards||[]).join("")}
  ${t.model?`<div class="foot">${esc(t.model)} · ${t.ms} ms · tool call → view${t.tests?` + ${t.tests} test${t.tests>1?"s":""}`:""}; numbers computed in your browser</div>`:""}</div>`).join("");
 const EXQ=["show NE cells near AR tumor","where do CHGA+ cells sit relative to AR-high regions?","compare SYP in core 3 vs core 7","do KLK3 and CHGA spots co-occur?"];
 const nokey=health&&!health.key_configured?`<div class="keynote"><b>Live LLM is off on this deployment.</b> Add <code>ANTHROPIC_API_KEY</code> in Vercel → project → Settings → Environment Variables, then redeploy. Until then questions use the in-browser rule parser, and the tests below still run on the data.</div>`:"";
 const direct=`<div class="dt"><div class="dth">Run a test directly</div><div class="chips">${DIRECT.map((d,i)=>`<button class="dtb" data-i="${i}">${esc(d[0])}</button>`).join("")}</div></div>`;
 const empty=turns.length?"":`<div class="empty"><div class="ic">✦</div><div>Ask anything about this section. A live LLM turns it into view changes and picks a test; the numbers are computed on the data in your browser.</div><div class="chips">${EXQ.map(e=>`<button>${e}</button>`).join("")}</div></div>`;
 return `<div class="llmbar">${st}<span class="hint">Enter in the top bar · shared links re-run through the LLM</span></div>${nokey}${empty}${T}${direct}`}
const DIRECT=[["NE-high vs AR-high: proximity",{test:"proximity",group_a:{kind:"program",value:"NE-high"},group_b:{kind:"program",value:"AR-high"}}],
 ["NE-high vs AR-high: neighbourhood enrichment",{test:"nhood_enrichment",group_a:{kind:"program",value:"NE-high"},group_b:{kind:"program",value:"AR-high"}}],
 ["CHGA > 1 near AR-high: co-occurrence",{test:"co_occurrence",group_a:{kind:"gene",gene:"CHGA",op:">",threshold:1},group_b:{kind:"program",value:"AR-high"},distance_um:300}],
 ["SYP, Leiden 3 vs 7: Wilcoxon",{test:"wilcoxon",group_a:{kind:"leiden",value:"3"},group_b:{kind:"leiden",value:"7"},gene:"SYP"}]];
async function directTest(i){const [label,spec]=DIRECT[i];const t={q:label,status:"done",answer:"Direct test, no LLM involved: the groups and parameters are the preset ones shown on the card.",cards:[pendingCard(spec)],tests:1};turns.push(t);
 t.notes=await applyView(spec.test==="wilcoxon"?{color_by:{kind:"gene",gene:spec.gene}}:{zoom_to:"all"},[spec]);paint(label);
 try{const rec=await runTest(spec);t.cards[0]=card(rec)}catch(e){t.cards[0]=`<div class="rc"><div class="err">${esc(e.message)}</div></div>`}paint(label);window.__done=(window.__done||0)+1}
function paint(q){A.setAsk(q,panel());const el=document.querySelectorAll("#side .turn");if(el.length)el[el.length-1].scrollIntoView({block:"start"})}
async function llmAsk(q){q=(q||"").trim();if(!q)return;if(health&&!health.key_configured){paint(q);return A.ruleAsk(q)}
 const t={q,status:"thinking",cards:[]};turns.push(t);paint(q);
 const history=turns.slice(-4,-1).flatMap(x=>[{role:"user",content:x.q},...(x.answer?[{role:"assistant",content:x.answer}]:[])]);
 let j;try{const r=await fetch(ENDPOINT,{method:"POST",headers:{"content-type":"application/json"},body:JSON.stringify({q,history,view:{preset:A.state.preset,color:A.state.enc=="geneSelection"?A.state.gene:"sets"}})});
  j=await r.json().catch(()=>({error:"bad_json"}));
  if(!r.ok){t.status="error";if(j.error==="no_key"){health={key_configured:false};turns.pop();paint(q);return A.ruleAsk(q)}
   t.err=j.error==="rate_limited"?`Rate limit reached; try again in ${j.retry_after_s}s.`:`LLM call failed (${j.error||r.status}). Falling back to the in-browser rule parser.`;paint(q);if(j.error!=="rate_limited")setTimeout(()=>A.ruleAsk(q),50);return}}
 catch(e){t.status="error";t.err="Could not reach the LLM endpoint. Falling back to the in-browser rule parser.";paint(q);setTimeout(()=>A.ruleAsk(q),50);return}
 t.status="done";t.answer=j.answer;t.interp=j.interpretation;t.model=j.model;t.ms=j.ms;const tests=(j.tests||[]).slice(0,3);t.tests=tests.length;window.__lastLLM=j;
 t.notes=await applyView(j.view,tests);t.cards=tests.map(pendingCard);paint(q);
 for(const [i,x] of tests.entries()){try{const rec=await runTest(x);rec.turn=turns.length;t.cards[i]=card(rec)}catch(e){t.cards[i]=`<div class="rc"><div class="rch"><span>${esc(x.test)}</span><span class="tag">NOT RUN</span></div><div class="err">${esc(e.message)}</div></div>`}paint(q)}
 window.__done=(window.__done||0)+1}

/* distance scale: calibrated to the Visium array itself (spot centres are 100 µm apart), from the median nearest-neighbour spacing in image pixels */
export function pitchPx(x,y){const n=x.length,d=new Float64Array(n);for(let i=0;i<n;i++){let b=Infinity;for(let j=0;j<n;j++){if(i===j)continue;const dx=x[i]-x[j],dy=y[i]-y[j],q=dx*dx+dy*dy;if(q<b)b=q}d[i]=Math.sqrt(b)}d.sort();return n%2?d[(n-1)/2]:(d[n/2-1]+d[n/2])/2}
let UM=null,DAR=null;
export async function init(app){A=app;N=A.SP.barcode.length;const pp=pitchPx(A.MT.x,A.MT.y);UM=100/pp;XS=Float64Array.from(A.MT.x,v=>v*UM);YS=Float64Array.from(A.MT.y,v=>v*UM);
 window.__scale={pitch_px:pp,um_per_px:UM,v1_um_per_px:A.MT.um_per_px,v1_pitch_um:pp*A.MT.um_per_px};
 worker=new Worker(new URL("./worker.js",import.meta.url),{type:"module"});worker.onmessage=e=>{const p=pending.get(e.data.id);if(!p)return;pending.delete(e.data.id);e.data.ok?p.res(e.data.out):p.rej(new Error(e.data.error))};
 document.getElementById("side").addEventListener("click",e=>{const d=e.target.closest(".dtb");if(d){directTest(+d.dataset.i);return}const b=e.target.closest(".rcx");if(!b)return;const r=window.__results.find(x=>x.id===b.dataset.id);if(!r)return;const o={...r,result:{...r.result,null:undefined,dists:undefined},idxA:undefined,idxB:undefined};
  navigator.clipboard.writeText(JSON.stringify(o,null,1)).then(()=>A.toast("Result JSON copied"),()=>{})});
 window.__health=fetch(ENDPOINT,{method:"GET"}).then(r=>r.json()).then(h=>{health=h}).catch(()=>{health={key_configured:false,unreachable:true}}).finally(()=>{if(!turns.length&&A.isAskIdle())paint("")});
 window.__viewState=()=>{const S=A.state;return {preset:S.preset,enc:S.enc,gene:S.gene,sel:S.sel,extra:S.extra?S.extra.tree[0].children.map(c=>c.name+":"+c.set.length):null,view:S.view}};
 window.llmAsk=llmAsk;window.runTestDirect=async t=>card(await runTest(t));return {ask:llmAsk,panel}}
