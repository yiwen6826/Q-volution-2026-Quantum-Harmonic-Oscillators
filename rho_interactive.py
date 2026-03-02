import numpy as np
import plotly.graph_objects as go

# ==============================
# PARAMETERS
# ==============================

K_vals = np.linspace(-5.0, 5.0, 11)
B_vals = np.linspace(-5.0, 5.0, 21)

NA = len(A_vals)
NB = len(B_vals)

# ==============================
# PRECOMPUTE DATA
# ==============================

def f(a, b):
    x = np.linspace(-5, 5, 200)
    y = a * x**2 + b * x
    return x, y

traces = []
for i, a in enumerate(A_vals):
    for j, b in enumerate(B_vals):
        x, y = f(a, b)
        traces.append(
            go.Scatter(
                x=x,
                y=y,
                visible=False,
                name=f"a={a:.2f}, b={b:.2f}",
            )
        )

# Initial indices
i0 = NA // 2
j0 = NB // 2
traces[i0 * NB + j0].visible = True

fig = go.Figure(traces)
fig.update_layout(
    title="f(x) = a x² + b x",
    xaxis_title="x",
    yaxis_title="y",
)

# ==============================
# SLIDERS (DUMB — JS WILL FIX LOGIC)
# ==============================

slider_a = {
    "active": i0,
    "y": -0.15,
    "x": 0.1,
    "len": 0.8,
    "currentvalue": {"prefix": "a = "},
    "steps": [
        {
            "label": f"{a:.2f}",
            "method": "skip",  # handled by JS
        }
        for a in A_vals
    ],
}

slider_b = {
    "active": j0,
    "y": -0.30,
    "x": 0.1,
    "len": 0.8,
    "currentvalue": {"prefix": "b = "},
    "steps": [
        {
            "label": f"{b:.2f}",
            "method": "skip",
        }
        for b in B_vals
    ],
}

fig.update_layout(sliders=[slider_a, slider_b])

# ==============================
# EXPORT HTML
# ==============================

html = fig.to_html(
    include_plotlyjs="cdn",
    full_html=True,
    div_id="plot",
)

# ==============================
# INJECT JAVASCRIPT
# ==============================

js = f"""
<script>
const NA = {NA};
const NB = {NB};

let a_idx = {i0};
let b_idx = {j0};

const plot = document.getElementById("plot");

function update_visibility() {{
    const vis = Array(NA * NB).fill(false);
    const idx = a_idx * NB + b_idx;
    vis[idx] = true;
    Plotly.restyle(plot, {{visible: vis}});
}}

plot.on('plotly_sliderchange', function(e) {{
    const sliderIndex = e.slider._index;
    const stepIndex = e.step._index;

    if (sliderIndex === 0) {{
        a_idx = stepIndex;
    }} else if (sliderIndex === 1) {{
        b_idx = stepIndex;
    }}

    update_visibility();
}});
</script>
"""

html = html.replace("</body>", js + "\n</body>")

# ==============================
# CREATE HTML
# ==============================

with open("interactive_plot.html", "w") as f:
    f.write(html)

print("Saved interactive_plot.html")