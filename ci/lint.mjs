// Lint: HTML structure via htmlhint + JS syntax of every inline module/classic script via node --check.
import {readFileSync,writeFileSync,mkdtempSync} from "fs";import {execFileSync} from "child_process";import {join} from "path";import {tmpdir} from "os";
const pages=["index.html","atlas.html","showcase.html","brief.html","explorer.html"];let bad=0;const dir=mkdtempSync(join(tmpdir(),"lint-"));
for(const p of pages){const h=readFileSync(p,"utf8");let i=0;
 for(const m of h.matchAll(/<script(?![^>]*\bsrc=)([^>]*)>([\s\S]*?)<\/script>/g)){if(/importmap|application\/json/.test(m[1]))continue;
  const f=join(dir,`${p}.${i++}.${/module/.test(m[1])?"mjs":"cjs"}`);writeFileSync(f,m[2]);
  try{execFileSync(process.execPath,["--check",f],{stdio:"pipe"})}catch(e){bad++;console.error(`✗ ${p} script #${i}: ${String(e.stderr).split("\n").slice(0,5).join("\n")}`)}}
 console.log(`✓ ${p}: ${i} inline script(s) parsed`)}
try{execFileSync("npx",["--yes","htmlhint@1.1.4","--config","ci/htmlhint.json",...pages],{stdio:"inherit"})}catch{bad++}
if(bad){console.error(`${bad} lint failure(s)`);process.exit(1)}
