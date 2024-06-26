#!/usr/bin/env python3

import json
import os
import base64

# Assuming these variables are set properly in your environment or script
png = "$png"
title = "$title"
description = "$description"

# Ensure the png path is an absolute path
absolute_path = os.path.abspath(png)

filename = os.path.basename(absolute_path)
name = filename.rstrip("_mcq.png")

id = title.lower().replace(" ", "_")

# Convert the PNG file to a base64 string
with open(absolute_path, "rb") as image_file:
    base64_string = base64.b64encode(image_file.read()).decode('utf-8')

# Create the data structure for JSON
data = {
    "id": f"benchmarking_{name}",
    "parent_id": id,
    "parent_name": title,
    "parent_description": description,
    "section_name": name,
    "plot_type": "image",
    "data": f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{base64_string}"/></div>'
}

# Write the JSON data to a file
with open(f"{name}_mqc.json", "w") as f:
    json.dump(data, f, indent=4)

print(f"JSON file created for {absolute_path}")

