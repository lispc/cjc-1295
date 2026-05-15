#!/usr/bin/env python3
"""Combine PyMOL panels with the data figure into a publication-ready composite."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# Load images
dala = Image.open("/home/scroll/personal/cjc-1295/workspace/step3/panel_DALA.png")
lala = Image.open("/home/scroll/personal/cjc-1295/workspace/step3/panel_LALA.png")
data_fig = Image.open("/home/scroll/personal/cjc-1295/workspace/results/mechanism_figure.png")

# Get dimensions
dw, dh = dala.size   # PyMOL panels
lw, lh = lala.size
fw, fh = data_fig.size  # mechanism figure (data + schematic)

# Resize PyMOL panels to same height as mechanism figure
target_h = fh
new_w_dala = int(dw * target_h / dh)
new_w_lala = int(lw * target_h / lh)
dala_rs = dala.resize((new_w_dala, target_h), Image.LANCZOS)
lala_rs = lala.resize((new_w_lala, target_h), Image.LANCZOS)

# Create composite: [DALA panel] [LALA panel] [Data figure]
# Top row: two PyMOL panels side by side
# Bottom: mechanism data figure
top_w = new_w_dala + new_w_lala
total_w = max(top_w, fw)
total_h = target_h + fh + 60  # 60px for gap

composite = Image.new('RGBA', (total_w, total_h), (255, 255, 255, 255))

# Place PyMOL panels on top row (centered)
top_offset_x = (total_w - top_w) // 2
composite.paste(dala_rs, (top_offset_x, 40))
composite.paste(lala_rs, (top_offset_x + new_w_dala, 40))

# Place data figure below
data_offset_x = (total_w - fw) // 2
composite.paste(data_fig, (data_offset_x, target_h + 60))

# Add labels
draw = ImageDraw.Draw(composite)
try:
    font_big = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 32)
    font_med = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
except:
    font_big = ImageFont.load_default()
    font_med = ImageFont.load_default()

# Title
draw.text((total_w//2, 5), "D-Ala2 Inhibits DPP-IV by Locking a Non-Productive Pose",
          fill=(30, 30, 30), font=font_big, anchor="mt")

# Panel labels
draw.text((top_offset_x + 10, 10), "A. D-Ala2 (CJC-1295) — tight but wrong angle",
          fill=(200, 40, 40), font=font_med, anchor="lt")
draw.text((top_offset_x + new_w_dala + 10, 10), "B. L-Ala2 (native GHRH) — productive angle but dissociates",
          fill=(40, 70, 140), font=font_med, anchor="lt")

# Data panel label
draw.text((data_offset_x + 10, target_h + 50), "C. Catalytic Geometry & Energetics",
          fill=(30, 30, 30), font=font_med, anchor="lt")

# Save
outpath = "/home/scroll/personal/cjc-1295/workspace/results/mechanism_composite.png"
composite.save(outpath, dpi=(300, 300))
print(f"Saved: {outpath}")
print(f"Dimensions: {composite.size}")
