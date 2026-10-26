// Vercel serverless function: POST /api/ask  ->  { answer, interpretation, view, tests, model }
// The model only chooses view changes and which test to run with which parameters.
// Every number shown in the app is computed afterwards by code in the browser (v2/stats.js), never by the model.
// The Anthropic key lives only in the Vercel env var ANTHROPIC_API_KEY. It is never sent to the client.

const MODELS = [process.env.ANTHROPIC_MODEL, "claude-sonnet-4-5", "claude-haiku-4-5", "claude-3-7-sonnet-latest"].filter(Boolean);
const API = (process.env.ANTHROPIC_BASE_URL || "https://api.anthropic.com").replace(/\/$/, "");
const LIMITS = { perMinute: +(process.env.RL_PER_MINUTE || 8), perDay: +(process.env.RL_PER_DAY || 120), globalPerHour: +(process.env.RL_GLOBAL_PER_HOUR || 600) };
const MAX_Q = 600;

// ---- rate limiting -------------------------------------------------------
// Upstash Redis (if UPSTASH_REDIS_REST_URL/TOKEN are set) gives a limit shared across instances.
// Otherwise an in-memory limiter applies per warm instance (best effort).
const mem = new Map();
function memHit(key, windowMs, limit) {
  const now = Date.now(); let e = mem.get(key);
  if (!e || now > e.reset) { e = { n: 0, reset: now + windowMs }; mem.set(key, e); }
  e.n++; if (mem.size > 5000) for (const [k, v] of mem) if (now > v.reset) mem.delete(k);
  return { ok: e.n <= limit, retry: Math.ceil((e.reset - now) / 1000) };
}
async function redisHit(key, windowS, limit) {
  const url = process.env.UPSTASH_REDIS_REST_URL, tok = process.env.UPSTASH_REDIS_REST_TOKEN;
  const r = await fetch(`${url}/pipeline`, { method: "POST", headers: { Authorization: `Bearer ${tok}`, "content-type": "application/json" },
    body: JSON.stringify([["INCR", key], ["EXPIRE", key, String(windowS), "NX"], ["TTL", key]]) });
  const j = await r.json(); const n = j[0].result, ttl = j[2].result;
  return { ok: n <= limit, retry: ttl > 0 ? ttl : windowS };
}
async function limit(ip) {
  const useRedis = !!(process.env.UPSTASH_REDIS_REST_URL && process.env.UPSTASH_REDIS_REST_TOKEN);
  const hit = useRedis ? (k, ms, l) => redisHit("rl:" + k, Math.ceil(ms / 1000), l) : async (k, ms, l) => memHit(k, ms, l);
  const checks = [await hit(`m:${ip}:${Math.floor(Date.now() / 60000)}`, 60000, LIMITS.perMinute),
                  await hit(`d:${ip}:${new Date().toISOString().slice(0, 10)}`, 86400000, LIMITS.perDay),
                  await hit(`g:${Math.floor(Date.now() / 3600000)}`, 3600000, LIMITS.globalPerHour)];
  const bad = checks.find(c => !c.ok);
  return { ok: !bad, retry: bad ? bad.retry : 0, store: useRedis ? "redis" : "memory" };
}

// ---- prompt + tool --------------------------------------------------------
const SYSTEM = `You drive an interactive spatial transcriptomics viewer (Vitessce) for ONE tissue section: GEO GSE230282 / GSM7211257, a 10x Visium CytAssist FFPE human prostate sample where de novo neuroendocrine prostate cancer (NE) sits next to AR-driven adenocarcinoma (Watanabe et al., IJMS 2023).

Facts you must respect:
- The unit is the Visium SPOT (55 um spot, 100 um centre-to-centre), not the single cell. 4,394 spots, 18,939 genes, log-normalised expression. If the user says "cells", read it as spots and say so in "interpretation".
- There is one section and one patient. There are no TMA cores. If the user says "core N", map it to Leiden cluster N (0-11) and say so in "interpretation".
- Spot groupings available (GroupSpec.kind):
  * "program": AR-high (628 spots, AR program z>1), NE-high (890, NE program z>1), "Dual AR+/NE+" (22), "Low/other" (2,854). Best default for "AR tumour" and "NE" questions.
  * "state": abstention-gated ensemble call: "AR adenocarcinoma", "NE / t-SCNC-like", "Amphicrine (AR+NE)", "Stroma", "Other", "Abstained".
  * "region": "AR region", "AR boundary", "Dual AR+/NE+", "Stroma/other", "NE boundary", "NE region".
  * "leiden": clusters "0" to "11".
  * "gene": spots where a gene passes a threshold, e.g. {"kind":"gene","gene":"CHGA","op":">","threshold":1}. Use HGNC symbols (chromogranin A = CHGA, synaptophysin = SYP, PSA = KLK3).
  * "score": {"kind":"score","score":"AR" or "NE","op":">","threshold":1} on the program z-scores.
  * "distance_to_ar": {"kind":"distance_to_ar","op":"<=","threshold":500} micrometres to the nearest AR-high spot.
  * "rest": every spot not in the other group. Any spec may add "and":[GroupSpec,...] to intersect.
- Tests you can request (the app computes them in code; you never compute or guess numbers):
  * "proximity": is group A closer to group B than label shuffling predicts? (mean nearest A-to-B distance, label-permutation null; null_scope "ab" shuffles only among A and B spots, a stricter null that keeps the tumour footprint fixed)
  * "nhood_enrichment": squidpy-style neighbourhood enrichment z-score of A-B contacts on a radius graph, label-permutation null.
  * "co_occurrence": squidpy co-occurrence ratio p(A | near B) / p(A) across distances, with a permutation p-value at distance_um.
  * "wilcoxon": differential expression of one gene between A and B (Wilcoxon rank-sum). Needs "gene".
- Presets: overview (cell states), arne (AR vs NE program scores), niche (AR/NE boundary), gene (one gene), conf (voter agreement).

Rules:
- Always call the respond tool exactly once.
- "view" is the NEW view for this question. The "Current view" line is context only: do not carry over a gene colouring unless this question is about that gene.
- When a test compares groups, highlight those two groups (A pink #ec4899, B blue #4f8cff) instead of colouring by a gene, except for wilcoxon, where colouring by the tested gene is better.
- "answer": 1-2 short plain sentences saying what you are showing and which test will run. Do NOT state or predict any statistic, p-value, direction or result. The app prints the computed numbers.
- No causal, mechanistic or clinical claims (no "drives", "causes", "predicts response", "prognosis"). Spot-level, one sample, descriptive.
- Prefer one test; use at most 3. Only request a test when the question asks something testable (near, closer, enriched, co-localise, compare, differ, higher in).
- For "compare GENE in X vs Y" use wilcoxon with that gene and also colour the view by the gene.
- For "where is GENE" colour by the gene; no test needed.
- Keep user words like "near" as a proximity or nhood test, and "co-localise / co-occur / sit together" as co_occurrence.
- No em dashes.`;

const GROUP = { type: "object", properties: {
  kind: { type: "string", enum: ["program", "state", "region", "leiden", "gene", "score", "distance_to_ar", "rest"] },
  value: { type: "string" }, gene: { type: "string" }, score: { type: "string", enum: ["AR", "NE"] },
  op: { type: "string", enum: [">", ">=", "<", "<="] }, threshold: { type: "number" },
  and: { type: "array", items: { type: "object" } }, label: { type: "string", description: "short display label" } }, required: ["kind"] };
const TOOL = { name: "respond", description: "Drive the viewer and choose statistical tests. Numbers are computed by the app, not by you.",
  input_schema: { type: "object", properties: {
    answer: { type: "string", description: "1-2 plain sentences, no result numbers" },
    interpretation: { type: "string", description: "how the question was mapped onto this data, e.g. cells read as spots" },
    view: { type: "object", properties: {
      preset: { type: "string", enum: ["overview", "arne", "niche", "gene", "conf"] },
      color_by: { type: "object", properties: { kind: { type: "string", enum: ["state", "region", "leiden", "gene", "score"] }, gene: { type: "string" }, score: { type: "string", enum: ["AR", "NE", "stroma", "confidence", "distance_to_ar"] } } },
      highlight: { type: "array", items: { type: "object", properties: { group: GROUP, color: { type: "string", description: "hex like #ec4899" } }, required: ["group"] } },
      zoom_to: { type: "string", description: "'all', 'a', 'b' (test groups), or 'highlight'" } } },
    tests: { type: "array", maxItems: 3, items: { type: "object", properties: {
      test: { type: "string", enum: ["proximity", "nhood_enrichment", "co_occurrence", "wilcoxon"] },
      group_a: GROUP, group_b: GROUP, gene: { type: "string" },
      n_perms: { type: "integer", minimum: 99, maximum: 2000 }, radius_um: { type: "number" }, distance_um: { type: "number" },
      null_scope: { type: "string", enum: ["all", "ab"], description: "proximity only: shuffle labels across all spots (default) or only among A and B spots" } },
      required: ["test", "group_a", "group_b"] } } },
    required: ["answer", "view", "tests"] } };

// ---- helpers --------------------------------------------------------------
function send(res, status, obj, extra = {}) {
  res.statusCode = status; res.setHeader("content-type", "application/json; charset=utf-8"); res.setHeader("cache-control", "no-store");
  for (const [k, v] of Object.entries(extra)) res.setHeader(k, v); res.end(JSON.stringify(obj));
}
function cors(req, res) {
  const origin = req.headers.origin; if (!origin) return true;
  const allowed = (process.env.ALLOWED_ORIGINS || "https://logan196.github.io").split(",").map(s => s.trim());
  let ok = false; try { const u = new URL(origin); ok = u.host === req.headers.host || allowed.includes(origin) || u.hostname === "localhost" || u.hostname === "127.0.0.1"; } catch {}
  if (ok) { res.setHeader("access-control-allow-origin", origin); res.setHeader("vary", "origin"); res.setHeader("access-control-allow-headers", "content-type"); res.setHeader("access-control-allow-methods", "POST, GET, OPTIONS"); }
  return ok;
}
async function readBody(req) {
  if (req.body && typeof req.body === "object") return req.body;
  if (typeof req.body === "string") return JSON.parse(req.body || "{}");
  const chunks = []; let size = 0; for await (const c of req) { size += c.length; if (size > 20000) throw new Error("body too large"); chunks.push(c); }
  return JSON.parse(Buffer.concat(chunks).toString("utf8") || "{}");
}
const clean = s => String(s || "").replace(/[\u0000-\u001f]/g, " ").slice(0, MAX_Q);

export default async function handler(req, res) {
  const okOrigin = cors(req, res);
  if (req.method === "OPTIONS") { res.statusCode = okOrigin ? 204 : 403; return res.end(); }
  if (req.method === "GET") return send(res, 200, { ok: true, key_configured: !!process.env.ANTHROPIC_API_KEY, model: MODELS[0], limits: LIMITS, store: process.env.UPSTASH_REDIS_REST_URL ? "redis" : "memory" });
  if (req.method !== "POST") return send(res, 405, { error: "method_not_allowed" });
  if (!okOrigin) return send(res, 403, { error: "origin_not_allowed" });
  if (!process.env.ANTHROPIC_API_KEY) return send(res, 503, { error: "no_key", message: "ANTHROPIC_API_KEY is not set on this Vercel project." });
  const ip = String(req.headers["x-forwarded-for"] || req.socket?.remoteAddress || "?").split(",")[0].trim();
  const rl = await limit(ip).catch(() => ({ ok: true, store: "error" }));
  if (!rl.ok) return send(res, 429, { error: "rate_limited", retry_after_s: rl.retry }, { "retry-after": String(rl.retry) });
  let body; try { body = await readBody(req); } catch { return send(res, 400, { error: "bad_body" }); }
  const q = clean(body.q); if (!q.trim()) return send(res, 400, { error: "empty_question" });
  const hist = (Array.isArray(body.history) ? body.history : []).slice(-6).map(m => ({ role: m.role === "assistant" ? "assistant" : "user", content: clean(m.content) })).filter(m => m.content);
  const ctx = clean(JSON.stringify(body.view || {})).slice(0, 300);
  const messages = [];
  for (const m of hist) { if (messages.length && messages[messages.length - 1].role === m.role) messages[messages.length - 1].content += "\n" + m.content; else messages.push({ ...m }); }
  if (messages.length && messages[0].role !== "user") messages.shift();
  const userMsg = `Current view: ${ctx}\nQuestion: ${q}`;
  if (messages.length && messages[messages.length - 1].role === "user") messages[messages.length - 1].content += "\n" + userMsg; else messages.push({ role: "user", content: userMsg });
  let lastErr = null;
  for (const model of MODELS) {
    const t0 = Date.now();
    try {
      const r = await fetch(`${API}/v1/messages`, { method: "POST", headers: { "x-api-key": process.env.ANTHROPIC_API_KEY, "anthropic-version": "2023-06-01", "content-type": "application/json" },
        body: JSON.stringify({ model, max_tokens: 900, temperature: 0, system: SYSTEM, tools: [TOOL], tool_choice: { type: "tool", name: "respond" }, messages }) });
      const j = await r.json().catch(() => ({}));
      if (r.status === 404 || (r.status === 400 && /model/i.test(JSON.stringify(j)))) { lastErr = { status: r.status, model }; continue; }
      if (!r.ok) return send(res, 502, { error: "llm_error", status: r.status, detail: (j.error && j.error.type) || "upstream" });
      const tu = (j.content || []).find(c => c.type === "tool_use" && c.name === "respond");
      if (!tu) return send(res, 502, { error: "no_tool_call" });
      return send(res, 200, { ...tu.input, model: j.model || model, ms: Date.now() - t0, usage: j.usage && { input_tokens: j.usage.input_tokens, output_tokens: j.usage.output_tokens }, rate: { store: rl.store } });
    } catch (e) { lastErr = { error: String(e.message || e).slice(0, 120), model }; }
  }
  return send(res, 502, { error: "llm_unavailable", detail: lastErr });
}
