"""Headless smoke test. Fails on any page error, console error or HTTP >=400
except the known zarr.js probes for `.zarray` inside zarr *groups* (harmless: zarr
checks for an array before falling back to the group's .zgroup).
Usage: python ci/smoke.py http://localhost:8765/  [--chrome /path/to/chrome]"""
import asyncio, re, sys, json
from playwright.async_api import async_playwright
BASE = sys.argv[1].rstrip("/") + "/"
CHROME = sys.argv[sys.argv.index("--chrome") + 1] if "--chrome" in sys.argv else None
ALLOW = re.compile(r"/vitessce/data\.zarr/(obs|var)(/[^/]+)?/\.zarray$")
ARGS = ["--no-sandbox", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]
CASES = [  # (path, viewport, check-js, wait ms)
    ("index.html?tour=0", (1440, 900), "window.__ready>0 && document.querySelectorAll('#vit canvas').length>=2 && document.querySelectorAll('#side .sec').length==4", 20000),
    ("index.html?preset=conf", (1440, 900), "document.querySelector('#presets .on').dataset.p=='conf'", 15000),
    ("index.html?q=compare%20KLK3%20in%20AR%20vs%20NE", (1440, 900), "/KLK3/.test(document.querySelector('#side .q')?.textContent) && /preset=arne/.test(location.search)", 20000),
    ("index.html?q=where%20are%20the%20immune%20cells", (1440, 900), "/PTPRC/.test(document.querySelector('#side').textContent)", 15000),
    ("index.html?preset=gene", (1440, 900), "/Distribution/.test(document.querySelector('#side').textContent)", 15000),
    ("index.html?q=where%20is%20chromogranin%20high", (1440, 900), "/CHGA/.test(document.querySelector('#side').textContent)", 20000),
    ("index.html", (390, 844), "document.documentElement.scrollWidth<=innerWidth+1 && document.querySelector('#boot').classList.contains('done')", 20000),
    ("atlas.html#/gene/KLK3", (1440, 900), "document.body.innerText.includes('KLK3')", 6000),
    ("showcase.html", (1440, 900), "document.querySelectorAll('canvas').length>0", 6000),
    ("brief.html", (1200, 900), "document.body.innerText.length>500", 2000),
    ("explorer.html", (1440, 900), "document.body.innerText.length>100", 5000),
]
async def main():
    fails, report = [], []
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path=CHROME, args=ARGS)
        for path, vp, check, wait in CASES:
            pg = await b.new_page(viewport={"width": vp[0], "height": vp[1]})
            errs = []; xreq = [0]
            pg.on("request", lambda r: xreq.__setitem__(0, xreq[0] + 1) if "/data.zarr/X/" in r.url else None)
            pg.on("pageerror", lambda e: errs.append("pageerror: " + str(e)[:300]))
            pg.on("console", lambda m: errs.append("console: " + m.text[:300]) if m.type == "error" and not m.text.startswith("Failed to load resource") else None)
            pg.on("response", lambda r: errs.append(f"http {r.status}: {r.url}") if r.status >= 400 and not ALLOW.search(r.url.split("?")[0]) else None)
            await pg.goto(BASE + path, wait_until="load", timeout=60000)
            await pg.wait_for_timeout(wait)
            ok = await pg.evaluate(f"()=>{{try{{return !!({check})}}catch(e){{return false}}}}")
            if not ok: errs.append("check failed: " + check)
            if xreq[0] > 40: errs.append(f"matrix budget: {xreq[0]} X-chunk requests (>40) - a view is loading the full matrix")
            report.append({"page": path, "viewport": vp, "ok": not errs, "errors": errs[:6]})
            if errs: fails.append(path)
            await pg.close()
        await b.close()
    print(json.dumps(report, indent=1))
    if fails: sys.exit(f"FAILED: {fails}")
    print("all smoke checks passed")
asyncio.run(main())
