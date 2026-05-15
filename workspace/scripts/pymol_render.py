#!/usr/bin/env python3
"""PyMOL rendering: D-Ala2 vs L-Ala2 catalytic geometry at DPP-IV active site."""
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

# Custom colors
cmd.set_color("dppiv_grey",  [0.72, 0.72, 0.75])
cmd.set_color("catalytic_c", [0.90, 0.45, 0.20])
cmd.set_color("oxyanion_c",  [0.18, 0.65, 0.55])
cmd.set_color("d_ala_c",     [0.90, 0.25, 0.30])
cmd.set_color("l_ala_c",     [0.28, 0.47, 0.62])
cmd.set_color("anchor_c",    [0.85, 0.30, 0.30])

# ── Load structures ──
cmd.load("frame_DALA_100ns.pdb", "DALA")
cmd.load("frame_LALA_25ns.pdb",  "LALA")

# Load starting structure as a 3rd reference
# cmd.load("frame_start.pdb", "start")

# ─────────────────────────────────────────────────────────
def setup_scene(obj_name, peptide_color, label_text):
    """Apply common styling to the complex."""
    cmd.disable("all")
    cmd.enable(obj_name)

    # DPP-IV: light cartoon
    cmd.show("cartoon", f"{obj_name} and chain A")
    cmd.color("dppiv_grey", f"{obj_name} and chain A")
    cmd.set("cartoon_transparency", 0.60, f"{obj_name} and chain A")

    # Catalytic triad Ser630, Asp708, His740
    cmd.show("sticks", f"{obj_name} and chain A and resi 630+708+740")
    cmd.color("catalytic_c", f"{obj_name} and chain A and resi 630+708+740")

    # Oxyanion hole: Tyr547, Ser631
    cmd.show("sticks", f"{obj_name} and chain A and resi 547+631")
    cmd.color("oxyanion_c", f"{obj_name} and chain A and resi 547+631")

    # Glu205/206 anchor
    cmd.show("sticks", f"{obj_name} and chain A and resi 205+206")
    cmd.color("anchor_c", f"{obj_name} and chain A and resi 205+206")

    # GHRH N-terminal peptide (resi 1-5): thick sticks
    cmd.show("sticks", f"{obj_name} and chain B and resi 1-5")
    cmd.color(peptide_color, f"{obj_name} and chain B and resi 1-5")
    cmd.set("stick_radius", 0.25, f"{obj_name} and chain B")

    # Highlight Ala2 CB (sidechain)
    cmd.show("spheres", f"{obj_name} and chain B and resi 2 and name CB")
    cmd.color(peptide_color, f"{obj_name} and chain B and resi 2 and name CB")
    cmd.set("sphere_scale", 0.45, f"{obj_name} and chain B and resi 2 and name CB")

    # Distance: Ser630 OG → Ala2 C
    cmd.distance(f"{obj_name}_dist",
                 f"{obj_name} and chain A and resi 630 and name OG",
                 f"{obj_name} and chain B and resi 2 and name C")
    cmd.color(peptide_color, f"{obj_name}_dist")
    cmd.set("dash_length", 0.15, f"{obj_name}_dist")
    cmd.set("dash_gap", 0.10, f"{obj_name}_dist")

    # Angle: OG – C – N
    cmd.angle(f"{obj_name}_ang",
              f"{obj_name} and chain A and resi 630 and name OG",
              f"{obj_name} and chain B and resi 2 and name C",
              f"{obj_name} and chain B and resi 2 and name N")
    cmd.color(peptide_color, f"{obj_name}_ang")

    # Labels
    cmd.set("label_size", 24)
    cmd.set("label_font_id", 5)
    cmd.set("label_outline_color", "white")

    # Zoom on active site
    cmd.zoom(f"{obj_name} and chain A and resi 630", 10)

    # Title label
    cmd.pseudoatom(f"{obj_name}_title", pos=[0, 0, 0])
    # Will be positioned after zoom


# ── Render D-Ala2 panel ──
setup_scene("DALA", "d_ala_c", "D-Ala2 (CJC-1295)")
# Reposition
cmd.turn("y", -20)
cmd.turn("x", 15)
cmd.move("z", 2)
cmd.ray(1600, 900)
cmd.png("panel_DALA.png", width=1600, height=900, dpi=300)
print("Panel DALA rendered")

# ── Render L-Ala2 panel ──
setup_scene("LALA", "l_ala_c", "L-Ala2 (native GHRH)")
cmd.turn("y", -20)
cmd.turn("x", 15)
cmd.move("z", 2)
cmd.ray(1600, 900)
cmd.png("panel_LALA.png", width=1600, height=900, dpi=300)
print("Panel LALA rendered")

print("Done — panel_DALA.png and panel_LALA.png")
