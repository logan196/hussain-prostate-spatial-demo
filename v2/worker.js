import {radiusGraph,nhoodEnrichment,proximity,coOccurrenceTest,wilcoxon} from "./stats.js";
/* runs the permutation tests off the main thread so the viewer stays responsive */
let G=null,GR=null;
self.onmessage=e=>{const m=e.data,t0=performance.now();
 try{let out;
  if(m.test==="nhood_enrichment"){if(!G||GR!==m.radius){G=radiusGraph(m.xs,m.ys,m.radius);GR=m.radius}
   out=nhoodEnrichment(G,m.labels,m.k,{nPerms:m.nPerms,seed:m.seed,pair:[0,1]});out.edges=G.edges;out.meanDegree=G.edges/G.n}
  else if(m.test==="proximity")out=proximity(m.xs,m.ys,m.labels,0,1,{nPerms:m.nPerms,seed:m.seed,withinUm:m.withinUm,scope:m.scope});
  else if(m.test==="co_occurrence")out=coOccurrenceTest(m.xs,m.ys,m.labels,m.k,0,1,m.thresholds,{distanceUm:m.distanceUm,nPerms:m.nPerms,seed:m.seed});
  else if(m.test==="wilcoxon")out=wilcoxon(m.a,m.b);
  out.ms=Math.round(performance.now()-t0);self.postMessage({id:m.id,ok:true,out})}
 catch(err){self.postMessage({id:m.id,ok:false,error:String(err&&err.message||err)})}};
