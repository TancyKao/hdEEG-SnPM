#!/usr/bin/env python3
"""Make a spectral/GLM HTML report self-contained for sharing.

The reports written by export_report.m / export_glm_report.m reference their
topoplot & periodogram PNGs by relative filename (the JS builds <img src="...">
at render time). Sharing the .html alone therefore shows broken images.

This bundles every PNG in the report folder into the HTML as base64 data URIs,
producing ONE file you can email / drop in Dropbox that renders anywhere offline.

Usage:
    python3 scripts/embed_report_images.py <report.html> [out.html]

Default out = <report>_standalone.html in the same folder. Stdlib only.
"""
import base64
import os
import sys
import json


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: embed_report_images.py <report.html> [out.html]")
    html_in = os.path.abspath(sys.argv[1])
    folder = os.path.dirname(html_in)
    out = sys.argv[2] if len(sys.argv) > 2 else \
        os.path.join(folder, os.path.splitext(os.path.basename(html_in))[0] + "_standalone.html")

    with open(html_in, "r", encoding="utf-8") as f:
        html = f.read()

    # base64-encode every PNG sitting next to the report, keyed by filename
    imgdata = {}
    total = 0
    for name in sorted(os.listdir(folder)):
        if name.lower().endswith(".png"):
            with open(os.path.join(folder, name), "rb") as fp:
                raw = fp.read()
            imgdata[name] = "data:image/png;base64," + base64.b64encode(raw).decode("ascii")
            total += len(raw)
    if not imgdata:
        sys.exit("No PNGs found next to %s" % html_in)

    # Inject the lookup table and make the <img> renderer prefer it.
    # The report builds images via:  <img src="${src}" ...>
    blob = "<script>window.IMGDATA=" + json.dumps(imgdata) + ";</script>\n"
    if "<img src=\"${src}\"" in html:
        html = html.replace('<img src="${src}"',
                             '<img src="${(window.IMGDATA&&IMGDATA[src])||src}"')
    else:
        print("warning: expected '<img src=\"${src}\"' not found; "
              "embedding table anyway (renderer may differ).", file=sys.stderr)
    # place the table just after <body ...> so it exists before any render
    lo = html.lower()
    i = lo.find("<body")
    if i != -1:
        j = html.find(">", i) + 1
        html = html[:j] + "\n" + blob + html[j:]
    else:
        html = blob + html

    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print("Embedded %d PNGs (%.1f MB raw) -> %s (%.1f MB)" %
          (len(imgdata), total / 1e6, out, os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main()
