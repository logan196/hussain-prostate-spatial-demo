import test from "node:test"; import assert from "node:assert/strict";
import {mulberry32,shuffle,radiusGraph,nhoodEnrichment,proximity,coOccurrence,wilcoxon,normSf} from "../v2/stats.js";
test("mulberry32 is deterministic", ()=>{const a=mulberry32(0),b=mulberry32(0);for(let i=0;i<5;i++)assert.equal(a(),b())});
test("normal tail matches reference values", ()=>{assert.ok(Math.abs(normSf(1.96).p-0.024997895148220435)<1e-12);assert.ok(Math.abs(normSf(10).log10p+23.118053405486076)<1e-6)});
test("wilcoxon matches scipy.stats.mannwhitneyu", ()=>{const r=wilcoxon([1,2,2,3,5,0,0],[0,0,1,1,2,0.5]);assert.equal(r.U1,29);assert.ok(Math.abs(r.p-0.2717921578243575)<1e-12)});
test("hex lattice: radius graph degree 6 inside, nhood z for a split lattice is strongly positive", ()=>{const xs=[],ys=[];for(let r=0;r<20;r++)for(let c=0;c<20;c++){xs.push(c*100+(r%2)*50);ys.push(r*86.6025)}
 const g=radiusGraph(xs,ys,130);assert.equal(g.indptr[10*20+10+1]-g.indptr[10*20+10],6);const lab=Int32Array.from(xs,x=>x<1000?0:1);const r=nhoodEnrichment(g,lab,2,{nPerms:200,seed:1,pair:[0,0]});assert.ok(r.zPair>10)});
test("proximity: two far clusters are farther than shuffled", ()=>{const xs=[],ys=[],lab=[];for(let i=0;i<200;i++){xs.push((i%20)*100);ys.push(Math.floor(i/20)*100);lab.push(i%20<3?0:i%20>16?1:2)}
 const r=proximity(xs,ys,Int32Array.from(lab),0,1,{nPerms:99,seed:0});assert.ok(r.pFarther<=0.02&&r.obs>1000)});
test("co-occurrence ratio is 1 for a single label", ()=>{const xs=[0,100,200,300],ys=[0,0,0,0];const o=coOccurrence(xs,ys,Int32Array.from([0,0,0,0]),1,[150,350]);assert.equal(o.at(0,0,0),1)});
