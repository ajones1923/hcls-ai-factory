#!/usr/bin/env python3
"""Generate the HCLS AI Factory architecture figure (single-glance SVG).

Color encodes the taxonomy — teal = compute engines, indigo = reasoning agents,
ember = the TSC disease-program vertical. Reproducible: edit the rosters below
and re-run. Output: docs/brief/architecture.svg
"""
from html import escape

# ---- palette (matches the capability brief) ----
INK        = "#171A1F"
INK_SOFT   = "#545C66"
GROUND     = "#F1F4F3"
PANEL      = "#FBFCFB"
LINE       = "#C6CECB"
ENGINE     = "#0C6E63"; ENGINE_SOFT = "#E4EFEC"; ENGINE_LINE = "#B8D6CF"
AGENT      = "#4C3F98"; AGENT_SOFT  = "#EBE9F5"; AGENT_LINE  = "#CFC9EA"
PROGRAM    = "#A8482B"; PROGRAM_SOFT= "#F5E7E0"; PROGRAM_LINE= "#E3C5B7"
PLATFORM   = "#232A31"

SERIF = "Charter, 'Iowan Old Style', Georgia, 'Times New Roman', serif"
SANS  = "system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
MONO  = "'SF Mono', 'JetBrains Mono', Menlo, Consolas, monospace"

engines = [
    ("E1", "Genomic Foundation", "5000"),
    ("E2", "Precision Intelligence", "5001"),
    ("E3", "Therapeutic Discovery", "8505"),
    ("E4", "Clinical Imaging", "8525"),
    ("E5", "Precision Oncology", "8526"),
    ("E6", "Cardiology", "8536"),
    ("E7", "Large-Molecule / Structural Biology", "8570–8578"),
    ("E8", "Single-Cell Analysis", "8573"),
]
agents = [
    ("CAR-T Intelligence", "8521"),
    ("Precision Biomarker Intelligence", "8528"),
    ("Pharmacogenomics Intelligence", "8507"),
    ("Precision Autoimmune Intelligence", "8531"),
    ("Neurology Intelligence", "8529"),
    ("Clinical Trial Intelligence", "8128"),
    ("Rare Disease Intelligence", "8544"),
    ("Single-Cell Intelligence", "8130"),
]
tsc_subagents = ["Variant Curator", "Trajectory Modeler", "Therapeutics Strategist",
                 "Phenome Mapper", "TAND Surveillance"]
platform = ["Capability Registry", "MCP Tool-Surface", "Workflow Composer",
            "Single-box MLOps", "Governance Gates"]

# ---- layout ----
W = 1440
MX = 44                      # outer margin x
INNER = W - 2 * MX
COLS = 4
GAP = 18
NODE_W = (INNER - (COLS - 1) * GAP) / COLS
NODE_H = 66
ROW_GAP = 16

parts = []

def wrap(text, maxchars):
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 <= maxchars or not cur:
            cur = (cur + " " + w).strip()
        else:
            lines.append(cur); cur = w
    if cur:
        lines.append(cur)
    return lines[:2]

def node(x, y, name, port, fill, stroke, accent, num=None):
    parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{NODE_W:.1f}" height="{NODE_H}" '
                 f'rx="5" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
    parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="4" height="{NODE_H}" rx="2" fill="{accent}"/>')
    tx = x + 16
    lines = wrap(name, 26)
    ty = y + (26 if len(lines) == 2 else 30)
    if num:
        parts.append(f'<text x="{tx}" y="{y+22}" font-family="{MONO}" font-size="12" '
                     f'font-weight="600" fill="{accent}">{num}</text>')
        tx_name = tx + 30
    else:
        tx_name = tx
    for i, ln in enumerate(lines):
        parts.append(f'<text x="{tx_name if i==0 else tx}" y="{ty + i*16}" font-family="{SANS}" '
                     f'font-size="14.5" font-weight="600" fill="{INK}">{escape(ln)}</text>')
    parts.append(f'<text x="{x + NODE_W - 14:.1f}" y="{y + NODE_H - 12}" text-anchor="end" '
                 f'font-family="{MONO}" font-size="11.5" fill="{INK_SOFT}">:{port}</text>')

def band_label(x, y, text, color):
    parts.append(f'<text x="{x}" y="{y}" font-family="{MONO}" font-size="12.5" '
                 f'letter-spacing="2.5" fill="{color}">{escape(text.upper())}</text>')

y = 0
# ===== Title =====
y = 58
parts.append(f'<text x="{MX}" y="{y}" font-family="{SERIF}" font-size="34" font-weight="600" '
             f'fill="{INK}">HCLS AI Factory</text>')
parts.append(f'<text x="{MX}" y="{y+30}" font-family="{SANS}" font-size="16" fill="{INK_SOFT}">'
             f'8 Engines &#183; 8 Intelligence Agents &#183; 1 Disease Program — patient DNA to therapeutic candidates, in hours, on one box.</text>')

# ===== Disease-program vertical (TSC) band =====
y += 66
band_h = 78
parts.append(f'<rect x="{MX}" y="{y}" width="{INNER}" height="{band_h}" rx="6" '
             f'fill="{PROGRAM_SOFT}" stroke="{PROGRAM_LINE}" stroke-width="1.5"/>')
parts.append(f'<rect x="{MX}" y="{y}" width="6" height="{band_h}" rx="3" fill="{PROGRAM}"/>')
band_label(MX+22, y+24, "Disease Program · the vertical", PROGRAM)
parts.append(f'<text x="{MX+22}" y="{y+50}" font-family="{SERIF}" font-size="20" font-weight="600" '
             f'fill="{INK}">Tuberous Sclerosis Complex</text>')
# sub-agents as pills on the right
sx = MX + 360
parts.append(f'<text x="{sx}" y="{y+30}" font-family="{MONO}" font-size="10.5" letter-spacing="1.5" '
             f'fill="{PROGRAM}">FIVE DISEASE-SPECIFIC AGENTS</text>')
px = sx
for sa in tsc_subagents:
    pw = 11 + len(sa) * 6.7
    parts.append(f'<rect x="{px:.1f}" y="{y+40}" width="{pw:.1f}" height="24" rx="12" '
                 f'fill="{PANEL}" stroke="{PROGRAM_LINE}"/>')
    parts.append(f'<text x="{px + pw/2:.1f}" y="{y+56}" text-anchor="middle" font-family="{SANS}" '
                 f'font-size="12" fill="{INK_SOFT}">{escape(sa)}</text>')
    px += pw + 8
parts.append(f'<text x="{W-MX}" y="{y+50}" text-anchor="end" font-family="{SERIF}" font-style="italic" '
             f'font-size="13" fill="{PROGRAM}">composes the horizontal layers ▼</text>')

# ===== Agents (reasoning) =====
y += band_h + 30
band_label(MX, y, "Intelligence Agents · the reasoning layer", AGENT)
y += 16
for i, (name, port) in enumerate(agents):
    col, row = i % COLS, i // COLS
    x = MX + col * (NODE_W + GAP)
    ny = y + row * (NODE_H + ROW_GAP)
    node(x, ny, name, port, AGENT_SOFT, AGENT_LINE, AGENT)
y = y + 2 * NODE_H + ROW_GAP

# ===== Engines (compute) =====
y += 34
band_label(MX, y, "Engines · the compute layer", ENGINE)
y += 16
for i, (num, name, port) in enumerate(engines):
    col, row = i % COLS, i // COLS
    x = MX + col * (NODE_W + GAP)
    ny = y + row * (NODE_H + ROW_GAP)
    node(x, ny, name, port, ENGINE_SOFT, ENGINE_LINE, ENGINE, num=num)
y = y + 2 * NODE_H + ROW_GAP

# ===== Platform foundation =====
y += 34
plat_h = 62
parts.append(f'<rect x="{MX}" y="{y}" width="{INNER}" height="{plat_h}" rx="6" fill="{PLATFORM}"/>')
parts.append(f'<text x="{MX+22}" y="{y+26}" font-family="{MONO}" font-size="11.5" letter-spacing="2" '
             f'fill="#9BB0AD">ONE PLATFORM · THE GLUE</text>')
px = MX + 22
for j, item in enumerate(platform):
    if j:
        parts.append(f'<text x="{px:.1f}" y="{y+48}" font-family="{SANS}" font-size="14" fill="#5C6B69">·</text>')
        px += 16
    parts.append(f'<text x="{px:.1f}" y="{y+48}" font-family="{SANS}" font-size="14.5" '
                 f'fill="#EAF0EF">{escape(item)}</text>')
    px += 9 + len(item) * 8.1
parts.append(f'<text x="{W-MX}" y="{y+40}" text-anchor="end" font-family="{SANS}" font-size="12.5" '
             f'fill="#8AA09D">one box · elastic burst to remote GPU</text>')
y += plat_h

# ===== Footer: legend + honesty =====
y += 40
def swatch(x, color, label):
    parts.append(f'<rect x="{x}" y="{y-11}" width="14" height="14" rx="3" fill="{color}"/>')
    parts.append(f'<text x="{x+20}" y="{y}" font-family="{SANS}" font-size="13" fill="{INK_SOFT}">{escape(label)}</text>')
swatch(MX, ENGINE, "Compute engine")
swatch(MX+170, AGENT, "Reasoning agent")
swatch(MX+345, PROGRAM, "Disease program")
parts.append(f'<text x="{W-MX}" y="{y}" text-anchor="end" font-family="{SANS}" font-size="12.5" '
             f'fill="{PROGRAM}">All clinical output is decision support for a qualified clinician — never autonomous diagnosis.</text>')

H = int(y + 34)
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
       f'font-family="{SANS}" role="img" aria-label="HCLS AI Factory architecture: 8 engines, 8 intelligence agents, 1 disease program on one platform">\n'
       f'<rect width="{W}" height="{H}" fill="{GROUND}"/>\n' + "\n".join(parts) + "\n</svg>\n")

with open("docs/brief/architecture.svg", "w") as f:
    f.write(svg)
print(f"wrote docs/brief/architecture.svg  ({W}x{H})")
