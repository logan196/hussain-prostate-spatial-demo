// Local dev server: static files with HTTP Range support (the OME-TIFF needs it) + /api/ask routed to the Vercel handler.
// Usage: ANTHROPIC_API_KEY=... node scripts/dev-server.mjs [port]
import http from "node:http"; import fs from "node:fs"; import path from "node:path"; import { fileURLToPath } from "node:url";
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), ".."); const PORT = +(process.argv[2] || 8765);
const handler = (await import(path.join(ROOT, "api/ask.js"))).default;
const MIME = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".mjs": "text/javascript", ".json": "application/json", ".css": "text/css", ".png": "image/png", ".jpg": "image/jpeg", ".tif": "image/tiff", ".pdf": "application/pdf", ".svg": "image/svg+xml", ".gif": "image/gif", ".md": "text/plain" };
http.createServer(async (req, res) => {
  const u = new URL(req.url, "http://x");
  if (u.pathname === "/api/ask") return handler(req, res);
  let p = decodeURIComponent(u.pathname); if (p.endsWith("/")) p += "index.html";
  const f = path.join(ROOT, p); if (!f.startsWith(ROOT) || !fs.existsSync(f) || fs.statSync(f).isDirectory()) { res.statusCode = 404; return res.end("not found"); }
  const size = fs.statSync(f).size, type = MIME[path.extname(f)] || "application/octet-stream";
  res.setHeader("accept-ranges", "bytes"); res.setHeader("content-type", type); res.setHeader("access-control-allow-origin", "*");
  const m = /bytes=(\d*)-(\d*)/.exec(req.headers.range || "");
  if (m) { const s = m[1] ? +m[1] : size - +m[2], e = m[1] && m[2] ? Math.min(+m[2], size - 1) : size - 1;
    res.statusCode = 206; res.setHeader("content-range", `bytes ${s}-${e}/${size}`); res.setHeader("content-length", e - s + 1); return fs.createReadStream(f, { start: s, end: e }).pipe(res); }
  res.setHeader("content-length", size); fs.createReadStream(f).pipe(res);
}).listen(PORT, () => console.log(`http://localhost:${PORT}/`));
