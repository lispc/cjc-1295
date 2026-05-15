#!/usr/bin/env python3
"""PyMOL: D-Ala2 vs L-Ala2 at DPP-IV active site — clean, no solvent."""
import sys, os
os.chdir("/home/scroll/personal/cjc-1295/workspace/step3")
sys.path.insert(0, "/home/scroll/miniforge3/lib/python3.13/site-packages")
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-c', '-q'])

# ── Settings ──
cmd.bg_color("white")
cmd.set("ray_trace_mode", 1)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("depth_cue", 0)
cmd.set("ray_opaque_background", 0)
cmd.set("cartoon_loop_radius", 0.3)
cmd.set("cartoon_oval_length", 1.8)
cmd.set("cartoon_rect_length", 1.8)

# Custom colors
cmd.set_color("dppiv_grey",  [0.68, 0.68, 0.72])
cmd.set_color("catalytic_c", [1.00, 0.50, 0.20])
cmd.set_color("oxyanion_c",  [0.15, 0.65, 0.55])
cmd.set_color("d_ala_c",     [0.90, 0.22, 0.28])
cmd.set_color("l_ala_c",     [0.22, 0.45, 0.62])
cmd.set_color("anchor_c",    [0.85, 0.30, 0.30])
cmd.set_color("pocket_c",    [0.90, 0.88, 0.70])

# ── Load structures (merged topology: DPP-IV=resi 39-766, GHRH=resi 1-29) ──
cmd.load("frame_DALA_clean.pdb", "DALA")
cmd.load("frame_LALA_clean.pdb",  "LALA")

# ── Selection helpers ──
DPPIV_SEL = "resi 39-766"
GHRH_SEL  = "resi 1-29"
CATALYTIC = "resi 630+708+740"
OXYANION  = "resi 547+631"
ANCHOR    = "resi 205+206"
ALA2_CB   = "resi 2 and name CB"
SER630_OG = "resi 630 and name OG"
ALA2_C    = "resi 2 and name C"
ALA2_N    = "resi 2 and name N"

def setup_scene(obj, pep_color):
    cmd.disable("all")
    cmd.enable(obj)

    # ── DPP-IV cartoon ──
    cmd.show("cartoon", f"{obj} and {DPPIV_SEL}")
    cmd.color("dppiv_grey", f"{obj} and {DPPIV_SEL}")
    cmd.set("cartoon_transparency", 0.50, f"{obj}")

    # ── Active site pocket surface (transparent) ──
    cmd.select("pocket_sel", f"{obj} and {DPPIV_SEL} and byres ({CATALYTIC} around 8)")
    cmd.show("surface", "pocket_sel")
    cmd.color("pocket_c", "pocket_sel")
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.55, "pocket_sel")

    # ── Catalytic triad (Ser630, Asp708, His740): orange sticks ──
    cmd.show("sticks", f"{obj} and {CATALYTIC}")
    cmd.color("catalytic_c", f"{obj} and {CATALYTIC}")
    cmd.set("stick_radius", 0.22, f"{obj} and {CATALYTIC}")

    # ── Oxyanion hole (Tyr547, Ser631): teal sticks ──
    cmd.show("sticks", f"{obj} and {OXYANION}")
    cmd.color("oxyanion_c", f"{obj} and {OXYANION}")

    # ── Glu205/206 anchor: red sticks ──
    cmd.show("sticks", f"{obj} and {ANCHOR}")
    cmd.color("anchor_c", f"{obj} and {ANCHOR}")

    # ── GHRH N-terminus (resi 1-5): thick sticks ──
    cmd.show("sticks", f"{obj} and {GHRH_SEL} and resi 1-5")
    cmd.color(pep_color, f"{obj} and {GHRH_SEL} and resi 1-5")
    cmd.set("stick_radius", 0.28, f"{obj} and {GHRH_SEL}")

    # ── GHRH rest (resi 6-29): thin cartoon ──
    cmd.show("cartoon", f"{obj} and {GHRH_SEL} and resi 6-29")
    cmd.color(pep_color, f"{obj} and {GHRH_SEL} and resi 6-29")
    cmd.set("cartoon_transparency", 0.30, f"{obj} and {GHRH_SEL}")

    # ── Ala2 CB sphere ──
    cmd.show("spheres", f"{obj} and {ALA2_CB}")
    cmd.color(pep_color, f"{obj} and {ALA2_CB}")
    cmd.set("sphere_scale", 0.50, f"{obj} and {ALA2_CB}")

    # ── Distance: Ser630 OG → Ala2 C ──
    cmd.distance(f"{obj}_d", f"{obj} and {SER630_OG}", f"{obj} and {ALA2_C}")
    cmd.color(pep_color, f"{obj}_d")
    cmd.set("dash_length", 0.12, f"{obj}_d")
    cmd.set("dash_gap", 0.08, f"{obj}_d")
    cmd.set("dash_radius", 0.08, f"{obj}_d")

    # ── Angle: OG – C – N ──
    cmd.angle(f"{obj}_a", f"{obj} and {SER630_OG}",
              f"{obj} and {ALA2_C}", f"{obj} and {ALA2_N}")
    cmd.color(pep_color, f"{obj}_a")

    # ── Labels ──
    cmd.set("label_size", 26)
    cmd.set("label_font_id", 5)
    cmd.set("label_outline_color", "white")
    cmd.set("label_color", pep_color)

    # ── Zoom on active site ──
    cmd.zoom(f"{obj} and {CATALYTIC}", 13)

    # ── Clean up temporary selection ──
    cmd.delete("pocket_sel")


# ── Render D-Ala2 ──
setup_scene("DALA", "d_ala_c")
cmd.turn("y", -25)
cmd.turn("x", 10)
cmd.move("z", 3)
cmd.move("y", -1)
cmd.ray(1600, 1000)
cmd.png("panel_DALA.png", width=1600, height=1000, dpi=300)
print("Panel DALA done")

# ── Render L-Ala2 ──
setup_scene("LALA", "l_ala_c")
cmd.turn("y", -25)
cmd.turn("x", 10)
cmd.move("z", 3)
cmd.move("y", -1)
cmd.ray(1600, 1000)
cmd.png("panel_LALA.png", width=1600, height=1000, dpi=300)
print("Panel LALA done")

print("Done.")
