// Playwright end-to-end run: 5 questions through the live LLM, view-change checks, screenshots, video, browser results for verify_tests.py
// Usage: CHROME=/path/to/chrome node verify/e2e.mjs http://localhost:8765/ out_dir

import { chromium } from "playwright-core";
import fs from "node:fs"; import crypto from "node:crypto";
const BASE=process.argv[2], OUT=process.argv[3]; fs.mkdirSync(OUT+"/shots",{recursive:true}); fs.mkdirSync(OUT+"/video",{recursive:true});
const Q=["show NE cells near AR tumor","where do CHGA+ cells sit relative to AR-high regions?","compare SYP in core 3 vs core 7","are AR-high and NE-high spots neighbours more often than chance?","where is KLK3 high?"];
const browser=await chromium.launch({executablePath:process.env.CHROME||undefined,args:["--use-gl=angle","--use-angle=swiftshader","--enable-unsafe-swiftshader","--ignore-gpu-blocklist"]});
const ctx=await browser.newContext({viewport:{width:1440,height:900},deviceScaleFactor:1,recordVideo:{dir:OUT+"/video",size:{width:1440,height:900}}});
const page=await ctx.newPage();const errs=[];page.on("pageerror",e=>errs.push("pageerror: "+e.message));page.on("console",m=>{if(m.type()==="error")errs.push("console: "+m.text().slice(0,200))});
const sha=async()=>{const b=await page.locator("#vit").screenshot();return crypto.createHash("sha1").update(b).digest("hex").slice(0,12)};
await page.goto(BASE+"index.html?tour=0",{waitUntil:"domcontentloaded"});
await page.waitForFunction(()=>window.__ready&&window.__llmReady,null,{timeout:120000});await page.waitForTimeout(2500);
const health=await page.evaluate(()=>window.__health&&window.__health.then(()=>document.querySelector(".llmpill")?.textContent));
const log=[];
for(const [i,q] of Q.entries()){
 const before=await page.evaluate(()=>window.__viewState()),hb=await sha(),done0=await page.evaluate(()=>window.__done||0),nres0=await page.evaluate(()=>window.__results.length);
 await page.click("#q");await page.fill("#q","");await page.type("#q",q,{delay:28});await page.keyboard.press("Enter");const t0=Date.now();
 await page.waitForFunction(d=>(window.__done||0)>d||document.querySelector("#side .turn:last-of-type .err"),done0,{timeout:180000});
 const tLLM=Date.now()-t0;await page.waitForTimeout(3500);
 const after=await page.evaluate(()=>window.__viewState()),ha=await sha(),llm=await page.evaluate(()=>window.__lastLLM);
 const res=await page.evaluate(n=>window.__results.slice(n).map(r=>({id:r.id,test:r.test,a:r.a,b:r.b,gene:r.gene,params:r.params,obs:r.result.obs??r.result.U1,z:r.result.z??r.result.zPair,p:{closer:r.result.pCloser,farther:r.result.pFarther,enriched:r.result.pEnriched,depleted:r.result.pDepleted,higher:r.result.pHigher,lower:r.result.pLower,two_sided:r.result.p}})),nres0);
 await page.screenshot({path:`${OUT}/shots/q${i+1}.png`});
 log.push({q,ms_to_done:tLLM,model:llm&&llm.model,llm_ms:llm&&llm.ms,answer:llm&&llm.answer,interpretation:llm&&llm.interpretation,view_call:llm&&llm.view,tests_called:llm&&llm.tests,before,after,view_changed:JSON.stringify(before)!==JSON.stringify(after),pixels_before:hb,pixels_after:ha,pixels_changed:hb!==ha,results:res});
 console.log("Q"+(i+1),JSON.stringify({q,ms:tLLM,changed:JSON.stringify(before)!==JSON.stringify(after),px:hb!==ha,tests:res.map(r=>r.test+" "+r.a.label+" vs "+r.b.label+" obs="+r.obs+" z="+r.z)}));
}
fs.writeFileSync(OUT+"/results_browser.json",JSON.stringify(await page.evaluate(()=>window.__results.map(r=>({...r,result:Object.fromEntries(Object.entries(r.result).filter(([k])=>k!=="null"&&k!=="dists"))})))));
// shared link: q= runs through the LLM
const p2=await ctx.newPage();await p2.goto(BASE+"index.html?tour=0&preset=overview&q=show+NE+cells+near+AR+tumor",{waitUntil:"domcontentloaded"});
await p2.waitForFunction(()=>(window.__done||0)>=1,null,{timeout:180000});await p2.waitForTimeout(3000);
const share=await p2.evaluate(()=>({llm:!!window.__lastLLM,model:window.__lastLLM&&window.__lastLLM.model,tests:window.__results.map(r=>r.test+": "+r.a.label+" vs "+r.b.label),view:window.__viewState(),url:location.href}));
await p2.screenshot({path:`${OUT}/shots/share_link.png`});console.log("SHARE",JSON.stringify(share));
fs.writeFileSync(OUT+"/e2e_log.json",JSON.stringify({base:BASE,health,questions:log,share,errors:errs},null,1));
console.log("ERRORS",JSON.stringify(errs.filter(e=>!/zarray|404/.test(e)).slice(0,10)));
const vpath=await page.video().path();await p2.close();await ctx.close();await browser.close();console.log("VIDEO",vpath);
