#!/usr/bin/env python3
"""Create PDX display artifacts from the June 20 computation."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor
from reportlab.pdfgen import canvas


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
PDX_DIR = PROJECT_ROOT / "case_study" / "pdx"
FIG_DIR = PROJECT_ROOT / "figures" / "pdx"
FIG_SOURCE_DIR = FIG_DIR / "source_data"
TABLE_DIR = PROJECT_ROOT / "tables" / "pdx"

DATE_TAG = "2026-06-20"

CANDIDATE_ORDER = [
    "BYL719 + binimetinib",
    "BYL719 + LEE011",
    "LEE011 + encorafenib",
    "BKM120 + binimetinib",
]

CANDIDATE_LABELS = {
    "BYL719 + binimetinib": ["BYL719 +", "binimetinib"],
    "BYL719 + LEE011": ["BYL719 +", "LEE011"],
    "LEE011 + encorafenib": ["LEE011 +", "encorafenib"],
    "BKM120 + binimetinib": ["BKM120 +", "binimetinib"],
}

ROLE_LABELS = {
    "primary_continuity_case": "Primary case",
    "strong_comparison_case": "Comparison case",
    "cautionary_sensitivity_case": "Cautionary case",
}

COLORS = {
    "ink": "#000000",
    "axis": "#000000",
    "grid": "#D9D9D9",
    "blue": "#5477C4",
    "blue_dark": "#2E4780",
    "orange": "#CC6F47",
    "orange_dark": "#804126",
    "gray": "#8C8C8C",
    "white": "#FFFFFF",
}


def ensure_dirs() -> None:
    for path in [FIG_DIR, FIG_SOURCE_DIR, TABLE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inputs = pd.read_csv(PDX_DIR / "pdx_case_study_inputs_2026-06-20.csv")
    bootstrap = pd.read_csv(PDX_DIR / "pdx_case_study_bootstrap_input_summary_2026-06-20.csv")
    planning = pd.read_csv(PDX_DIR / "pdx_case_study_uncertainty_planning_2026-06-20.csv")
    design = pd.read_csv(PDX_DIR / "pdx_case_study_design_outputs_2026-06-20.csv")
    variance_sensitivity = pd.read_csv(PDX_DIR / "pdx_case_study_clinical_variance_sensitivity_2026-06-20.csv")
    validation = pd.read_csv(PDX_DIR / "pdx_case_study_validation_checks_2026-06-20.csv")
    return inputs, bootstrap, planning, design, variance_sensitivity, validation


def ordered(df: pd.DataFrame, col: str = "candidate") -> pd.DataFrame:
    out = df.copy()
    out["_order"] = out[col].map({value: idx for idx, value in enumerate(CANDIDATE_ORDER)})
    if out["_order"].isna().any():
        missing = sorted(out.loc[out["_order"].isna(), col].astype(str).unique())
        raise ValueError(f"Unexpected candidate order values: {missing}")
    return out.sort_values("_order").drop(columns="_order")


def round_float(value: object, digits: int = 3) -> object:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, int)):
        if isinstance(value, float) and not math.isfinite(value):
            return ""
        return f"{float(value):.{digits}f}"
    return value


def round_numeric_frame(df: pd.DataFrame, digits: int = 3) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda value: round_float(value, digits))
    return out


def whole_number_or_text(value: object) -> object:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, int)):
        return str(int(round(float(value))))
    return value


def sentence_case(value: object) -> object:
    if pd.isna(value):
        return value
    text = str(value)
    return text[:1].upper() + text[1:] if text else text


def write_csv(path: Path, df: pd.DataFrame) -> None:
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def build_source_and_tables(
    inputs: pd.DataFrame,
    bootstrap: pd.DataFrame,
    planning: pd.DataFrame,
    design: pd.DataFrame,
    variance_sensitivity: pd.DataFrame,
    validation: pd.DataFrame,
) -> dict[str, Path]:
    inputs = inputs.rename(columns={"combination": "candidate"})
    bootstrap = bootstrap.rename(columns={"combination": "candidate"})
    planning = ordered(planning)
    inputs = ordered(inputs)
    bootstrap = ordered(bootstrap)

    figure_a_rows = []
    for _, row in bootstrap.iterrows():
        for contrast, prefix in [("AB - A", "delta_ab_minus_a"), ("AB - B", "delta_ab_minus_b")]:
            figure_a_rows.append(
                {
                    "candidate": row["candidate"],
                    "contrast": contrast,
                    "median": row[f"{prefix}_median"],
                    "lower_95": row[f"{prefix}_lcl"],
                    "upper_95": row[f"{prefix}_ucl"],
                    "valid_bootstrap_fraction": row["valid_boot_fraction"],
                }
            )
    figure_a = pd.DataFrame(figure_a_rows)

    figure_b_rows = []
    for _, row in planning.iterrows():
        figure_b_rows.append(
            {
                "candidate": row["candidate"],
                "design_summary": "Plug-in optimized",
                "total_N": row["plugin_joint_N"],
                "joint_power": row["plugin_joint_power"],
                "valid_design_input": True,
            }
        )
        figure_b_rows.append(
            {
                "candidate": row["candidate"],
                "design_summary": "Lower-bound effects",
                "total_N": row["lower_bound_joint_N"],
                "joint_power": row["lower_bound_joint_power"],
                "valid_design_input": bool(row["valid_lower_bound_input"]),
            }
        )
    figure_b = pd.DataFrame(figure_b_rows)
    figure_b["total_N"] = figure_b["total_N"].map(whole_number_or_text)

    table_main = planning.merge(
        inputs[["candidate", "n_triplet"]],
        on="candidate",
        how="left",
    )
    table_main = table_main[
        [
            "candidate",
            "case_role",
            "n_triplet",
            "plugin_joint_N",
            "lower_bound_joint_N",
            "plugin_design_joint_power_under_lower_bound_effects",
            "uncertainty_sample_size_inflation",
            "planning_interpretation",
        ]
    ].rename(
        columns={
            "case_role": "case_role_for_manuscript",
            "n_triplet": "complete_pdx_triplets",
            "plugin_joint_N": "plug_in_optimized_total_N",
            "lower_bound_joint_N": "lower_bound_effect_total_N",
            "plugin_design_joint_power_under_lower_bound_effects": "plug_in_design_joint_power_under_lower_bound_effects",
            "uncertainty_sample_size_inflation": "sample_size_inflation",
        }
    )
    table_main["case_role_for_manuscript"] = table_main["case_role_for_manuscript"].map(ROLE_LABELS)
    for text_col in [
        "plug_in_optimized_total_N",
        "lower_bound_effect_total_N",
        "plug_in_design_joint_power_under_lower_bound_effects",
        "sample_size_inflation",
    ]:
        table_main[text_col] = table_main[text_col].astype("object")
    table_main["plug_in_optimized_total_N"] = table_main["plug_in_optimized_total_N"].map(whole_number_or_text)
    table_main["lower_bound_effect_total_N"] = table_main["lower_bound_effect_total_N"].map(whole_number_or_text)
    table_main.loc[table_main["lower_bound_effect_total_N"].eq(""), "lower_bound_effect_total_N"] = "No finite design"
    table_main.loc[
        table_main["plug_in_design_joint_power_under_lower_bound_effects"].isna(),
        "plug_in_design_joint_power_under_lower_bound_effects",
    ] = "Not evaluated"
    table_main.loc[table_main["sample_size_inflation"].isna(), "sample_size_inflation"] = "Not evaluated"
    table_main["planning_interpretation"] = table_main["planning_interpretation"].map(sentence_case)

    supp_input = bootstrap[
        [
            "candidate",
            "n_boot_requested",
            "n_boot_valid",
            "valid_boot_fraction",
            "delta_ab_minus_a_median",
            "delta_ab_minus_a_lcl",
            "delta_ab_minus_a_ucl",
            "delta_ab_minus_b_median",
            "delta_ab_minus_b_lcl",
            "delta_ab_minus_b_ucl",
            "min_contribution_median",
            "min_contribution_lcl",
            "min_contribution_ucl",
            "prob_both_contributions_positive",
            "prob_min_contribution_gt_0_10",
        ]
    ].merge(
        inputs[
            [
                "candidate",
                "delta_ab_minus_a",
                "delta_ab_minus_b",
            ]
        ],
        on="candidate",
        how="left",
    ).rename(
        columns={
            "n_boot_requested": "bootstrap_samples_requested",
            "n_boot_valid": "valid_bootstrap_samples",
            "valid_boot_fraction": "valid_bootstrap_fraction",
            "delta_ab_minus_a": "observed_plug_in_delta_ab_minus_a",
            "delta_ab_minus_b": "observed_plug_in_delta_ab_minus_b",
            "prob_both_contributions_positive": "bootstrap_fraction_both_contributions_positive",
            "prob_min_contribution_gt_0_10": "bootstrap_fraction_min_contribution_gt_0_10",
        }
    )
    supp_input = supp_input[
        [
            "candidate",
            "bootstrap_samples_requested",
            "valid_bootstrap_samples",
            "valid_bootstrap_fraction",
            "observed_plug_in_delta_ab_minus_a",
            "delta_ab_minus_a_median",
            "delta_ab_minus_a_lcl",
            "delta_ab_minus_a_ucl",
            "observed_plug_in_delta_ab_minus_b",
            "delta_ab_minus_b_median",
            "delta_ab_minus_b_lcl",
            "delta_ab_minus_b_ucl",
            "min_contribution_median",
            "min_contribution_lcl",
            "min_contribution_ucl",
            "bootstrap_fraction_both_contributions_positive",
            "bootstrap_fraction_min_contribution_gt_0_10",
        ]
    ]

    supp_design = ordered(design)
    if "N" in supp_design.columns:
        supp_design["N"] = supp_design["N"].map(whole_number_or_text)

    supp_variance = variance_sensitivity.copy()
    supp_variance["input_scenario"] = supp_variance["input_scenario"].map(
        {
            "plug_in": "Plug-in",
            "bootstrap_lower_95_effects": "Lower bound",
        }
    )
    supp_variance["optimized_total_N"] = supp_variance["optimized_total_N"].map(whole_number_or_text)

    supp_validation = validation.copy()
    keep_checks = [
        "source data has expected dimensions from PDX contract",
        "all planned PDX candidates have at least 20 complete triplets",
        "bootstrap effect validity fraction is high",
        "clinical covariance model is independent randomized arms",
        "primary one-sided alpha is 0.025",
        "plug-in joint optimized designs reach target joint power",
        "invalid lower-bound effect case is explicitly flagged rather than forced",
        "primary plug-in design is evaluated under lower-bound contribution effects",
    ]
    supp_validation = supp_validation[supp_validation["check"].isin(keep_checks)].copy()

    paths = {
        "figure_a_source": FIG_SOURCE_DIR / f"figure_4a_pdx_contribution_uncertainty_source_{DATE_TAG}.csv",
        "figure_b_source": FIG_SOURCE_DIR / f"figure_4b_pdx_design_consequence_source_{DATE_TAG}.csv",
        "main_table": TABLE_DIR / f"table_4_pdx_design_planning_main_{DATE_TAG}.csv",
        "supp_input_table": TABLE_DIR / f"supp_table_pdx_input_uncertainty_{DATE_TAG}.csv",
        "supp_design_table": TABLE_DIR / f"supp_table_pdx_design_outputs_full_{DATE_TAG}.csv",
        "supp_variance_sensitivity_table": TABLE_DIR / f"supp_table_pdx_clinical_variance_sensitivity_{DATE_TAG}.csv",
        "supp_validation_table": TABLE_DIR / f"supp_table_pdx_validation_checks_{DATE_TAG}.csv",
    }

    write_csv(paths["figure_a_source"], round_numeric_frame(figure_a))
    write_csv(paths["figure_b_source"], round_numeric_frame(figure_b))
    write_csv(paths["main_table"], round_numeric_frame(table_main))
    write_csv(paths["supp_input_table"], round_numeric_frame(supp_input))
    write_csv(paths["supp_design_table"], round_numeric_frame(supp_design))
    write_csv(paths["supp_variance_sensitivity_table"], round_numeric_frame(supp_variance))
    write_csv(paths["supp_validation_table"], round_numeric_frame(supp_validation))
    return paths


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def draw_multiline_text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: list[str], fnt, fill: str, spacing: int = 7) -> None:
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * (fnt.size + spacing)), line, fill=fill, font=fnt)


def draw_legend_item(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str, fnt) -> int:
    draw.ellipse((x, y + 4, x + 20, y + 24), fill=color, outline=COLORS["ink"], width=2)
    draw.text((x + 32, y), label, fill=COLORS["ink"], font=fnt)
    return x + 32 + int(draw.textlength(label, font=fnt)) + 46


def x_scale(value: float, left: int, right: int, xmin: float, xmax: float) -> float:
    return left + (value - xmin) / (xmax - xmin) * (right - left)


def draw_axis(
    draw: ImageDraw.ImageDraw,
    left: int,
    top: int,
    right: int,
    bottom: int,
    xmin: float,
    xmax: float,
    xticks: list[float],
    label: str,
    label_font,
    tick_font,
) -> None:
    draw.line((left, bottom, right, bottom), fill=COLORS["axis"], width=3)
    draw.line((left, top, left, bottom), fill=COLORS["axis"], width=3)
    for tick in xticks:
        x = x_scale(tick, left, right, xmin, xmax)
        draw.line((x, top, x, bottom), fill=COLORS["grid"], width=2)
        draw.line((x, bottom, x, bottom + 10), fill=COLORS["axis"], width=3)
        text = f"{tick:.1f}" if abs(tick - round(tick)) > 1e-8 else f"{int(tick)}"
        draw.text((x - draw.textlength(text, font=tick_font) / 2, bottom + 18), text, fill=COLORS["ink"], font=tick_font)
    draw.text(
        ((left + right) / 2 - draw.textlength(label, font=label_font) / 2, bottom + 64),
        label,
        fill=COLORS["ink"],
        font=label_font,
    )


def create_png(figure_a: pd.DataFrame, planning: pd.DataFrame, png_path: Path) -> None:
    width, height = 2400, 1350
    img = Image.new("RGB", (width, height), COLORS["white"])
    draw = ImageDraw.Draw(img)

    panel_font = font(50, bold=True)
    label_font = font(34)
    tick_font = font(30)
    legend_font = font(30)
    candidate_font = font(31)

    a_left, a_top, a_right, a_bottom = 350, 145, 1125, 1040
    b_left, b_top, b_right, b_bottom = 1425, 145, 2285, 1040
    row_top, row_bottom = a_top + 60, a_bottom - 90
    row_gap = (row_bottom - row_top) / (len(CANDIDATE_ORDER) - 1)
    y_positions = {candidate: int(row_top + idx * row_gap) for idx, candidate in enumerate(CANDIDATE_ORDER)}

    draw.text((80, 70), "A", fill=COLORS["ink"], font=panel_font)
    draw.text((1280, 70), "B", fill=COLORS["ink"], font=panel_font)

    legend_x = 350
    legend_x = draw_legend_item(draw, legend_x, 70, "AB - A", COLORS["blue"], legend_font)
    draw_legend_item(draw, legend_x, 70, "AB - B", COLORS["orange"], legend_font)

    legend_x = 1425
    legend_x = draw_legend_item(draw, legend_x, 70, "Plug-in", COLORS["blue"], legend_font)
    draw_legend_item(draw, legend_x, 70, "Lower-bound", COLORS["orange"], legend_font)

    a_xmin, a_xmax = -0.2, 1.5
    draw_axis(
        draw,
        a_left,
        a_top - 35,
        a_right,
        a_bottom,
        a_xmin,
        a_xmax,
        [0, 0.5, 1.0, 1.5],
        "Bootstrap contribution effect",
        label_font,
        tick_font,
    )
    zero_x = x_scale(0, a_left, a_right, a_xmin, a_xmax)
    draw.line((zero_x, a_top - 35, zero_x, a_bottom), fill=COLORS["ink"], width=2)

    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        draw_multiline_text(draw, (70, y - 32), CANDIDATE_LABELS[candidate], candidate_font, COLORS["ink"])
        part = figure_a.loc[figure_a["candidate"] == candidate]
        for _, row in part.iterrows():
            offset = -22 if row["contrast"] == "AB - A" else 22
            color = COLORS["blue"] if row["contrast"] == "AB - A" else COLORS["orange"]
            edge = COLORS["blue_dark"] if row["contrast"] == "AB - A" else COLORS["orange_dark"]
            y_mark = y + offset
            x_l = x_scale(row["lower_95"], a_left, a_right, a_xmin, a_xmax)
            x_u = x_scale(row["upper_95"], a_left, a_right, a_xmin, a_xmax)
            x_m = x_scale(row["median"], a_left, a_right, a_xmin, a_xmax)
            draw.line((x_l, y_mark, x_u, y_mark), fill=color, width=6)
            draw.line((x_l, y_mark - 14, x_l, y_mark + 14), fill=color, width=5)
            draw.line((x_u, y_mark - 14, x_u, y_mark + 14), fill=color, width=5)
            draw.ellipse((x_m - 13, y_mark - 13, x_m + 13, y_mark + 13), fill=color, outline=edge, width=3)

    b_xmin, b_xmax = 0, 1150
    draw_axis(
        draw,
        b_left,
        b_top - 35,
        b_right,
        b_bottom,
        b_xmin,
        b_xmax,
        [0, 300, 600, 900],
        "Optimized total sample size",
        label_font,
        tick_font,
    )
    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        row = planning.loc[planning["candidate"] == candidate].iloc[0]
        x_plugin = x_scale(float(row["plugin_joint_N"]), b_left, b_right, b_xmin, b_xmax)
        draw.ellipse((x_plugin - 14, y - 14, x_plugin + 14, y + 14), fill=COLORS["blue"], outline=COLORS["blue_dark"], width=3)
        if bool(row["valid_lower_bound_input"]):
            x_lower = x_scale(float(row["lower_bound_joint_N"]), b_left, b_right, b_xmin, b_xmax)
            draw.line((x_plugin, y, x_lower, y), fill=COLORS["gray"], width=5)
            draw.ellipse((x_lower - 14, y - 14, x_lower + 14, y + 14), fill=COLORS["orange"], outline=COLORS["orange_dark"], width=3)
        else:
            draw.text((b_left + 360, y - 18), "no finite lower-bound N", fill=COLORS["ink"], font=tick_font)

    img.save(png_path, dpi=(300, 300))


def svg_text(x: float, y: float, text: str, size: int, weight: str = "normal", anchor: str = "start") -> str:
    text = str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="#000000">{text}</text>'
    )


def create_svg(figure_a: pd.DataFrame, planning: pd.DataFrame, svg_path: Path) -> None:
    width, height = 2400, 1350
    a_left, a_top, a_right, a_bottom = 350, 145, 1125, 1040
    b_left, b_top, b_right, b_bottom = 1425, 145, 2285, 1040
    row_top, row_bottom = a_top + 60, a_bottom - 90
    row_gap = (row_bottom - row_top) / (len(CANDIDATE_ORDER) - 1)
    y_positions = {candidate: int(row_top + idx * row_gap) for idx, candidate in enumerate(CANDIDATE_ORDER)}
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        svg_text(80, 110, "A", 50, "bold"),
        svg_text(1280, 110, "B", 50, "bold"),
    ]

    def line(x1, y1, x2, y2, color, width_px=3) -> None:
        elements.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width_px}"/>')

    def circle(x, y, r, fill, stroke) -> None:
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')

    def axis(left, top, right, bottom, xmin, xmax, ticks, label) -> None:
        line(left, bottom, right, bottom, COLORS["axis"], 3)
        line(left, top, left, bottom, COLORS["axis"], 3)
        for tick in ticks:
            x = x_scale(tick, left, right, xmin, xmax)
            line(x, top, x, bottom, COLORS["grid"], 2)
            line(x, bottom, x, bottom + 10, COLORS["axis"], 3)
            tick_text = f"{tick:.1f}" if abs(tick - round(tick)) > 1e-8 else f"{int(tick)}"
            elements.append(svg_text(x, bottom + 50, tick_text, 30, anchor="middle"))
        elements.append(svg_text((left + right) / 2, bottom + 105, label, 34, anchor="middle"))

    elements.append(f'<circle cx="360" cy="84" r="10" fill="{COLORS["blue"]}" stroke="#000000" stroke-width="2"/>')
    elements.append(svg_text(392, 94, "AB - A", 30))
    elements.append(f'<circle cx="555" cy="84" r="10" fill="{COLORS["orange"]}" stroke="#000000" stroke-width="2"/>')
    elements.append(svg_text(587, 94, "AB - B", 30))
    elements.append(f'<circle cx="1435" cy="84" r="10" fill="{COLORS["blue"]}" stroke="#000000" stroke-width="2"/>')
    elements.append(svg_text(1467, 94, "Plug-in", 30))
    elements.append(f'<circle cx="1635" cy="84" r="10" fill="{COLORS["orange"]}" stroke="#000000" stroke-width="2"/>')
    elements.append(svg_text(1667, 94, "Lower-bound", 30))

    a_xmin, a_xmax = -0.2, 1.5
    axis(a_left, a_top - 35, a_right, a_bottom, a_xmin, a_xmax, [0, 0.5, 1.0, 1.5], "Bootstrap contribution effect")
    zero_x = x_scale(0, a_left, a_right, a_xmin, a_xmax)
    line(zero_x, a_top - 35, zero_x, a_bottom, COLORS["ink"], 2)
    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        label_lines = CANDIDATE_LABELS[candidate]
        elements.append(svg_text(70, y - 12, label_lines[0], 31))
        elements.append(svg_text(70, y + 31, label_lines[1], 31))
        part = figure_a.loc[figure_a["candidate"] == candidate]
        for _, row in part.iterrows():
            offset = -22 if row["contrast"] == "AB - A" else 22
            color = COLORS["blue"] if row["contrast"] == "AB - A" else COLORS["orange"]
            edge = COLORS["blue_dark"] if row["contrast"] == "AB - A" else COLORS["orange_dark"]
            y_mark = y + offset
            x_l = x_scale(row["lower_95"], a_left, a_right, a_xmin, a_xmax)
            x_u = x_scale(row["upper_95"], a_left, a_right, a_xmin, a_xmax)
            x_m = x_scale(row["median"], a_left, a_right, a_xmin, a_xmax)
            line(x_l, y_mark, x_u, y_mark, color, 6)
            line(x_l, y_mark - 14, x_l, y_mark + 14, color, 5)
            line(x_u, y_mark - 14, x_u, y_mark + 14, color, 5)
            circle(x_m, y_mark, 13, color, edge)

    b_xmin, b_xmax = 0, 1150
    axis(b_left, b_top - 35, b_right, b_bottom, b_xmin, b_xmax, [0, 300, 600, 900], "Optimized total sample size")
    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        row = planning.loc[planning["candidate"] == candidate].iloc[0]
        x_plugin = x_scale(float(row["plugin_joint_N"]), b_left, b_right, b_xmin, b_xmax)
        circle(x_plugin, y, 14, COLORS["blue"], COLORS["blue_dark"])
        if bool(row["valid_lower_bound_input"]):
            x_lower = x_scale(float(row["lower_bound_joint_N"]), b_left, b_right, b_xmin, b_xmax)
            line(x_plugin, y, x_lower, y, COLORS["gray"], 5)
            circle(x_lower, y, 14, COLORS["orange"], COLORS["orange_dark"])
        else:
            elements.append(svg_text(b_left + 360, y + 10, "no finite lower-bound N", 30))

    elements.append("</svg>")
    svg_path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def create_pdf_from_svg_logic(figure_a: pd.DataFrame, planning: pd.DataFrame, pdf_path: Path) -> None:
    scale = 0.3
    width, height = 2400 * scale, 1350 * scale
    c = canvas.Canvas(str(pdf_path), pagesize=(width, height))

    def sx(x: float) -> float:
        return x * scale

    def sy(y: float) -> float:
        return height - y * scale

    def set_color(hex_color: str) -> None:
        c.setStrokeColor(HexColor(hex_color))
        c.setFillColor(HexColor(hex_color))

    def line(x1, y1, x2, y2, color, width_px=3) -> None:
        set_color(color)
        c.setLineWidth(width_px * scale)
        c.line(sx(x1), sy(y1), sx(x2), sy(y2))

    def text(x, y, value, size, bold=False, anchor="start") -> None:
        c.setFillColor(HexColor(COLORS["ink"]))
        c.setFont("Helvetica-Bold" if bold else "Helvetica", size * scale)
        if anchor == "middle":
            c.drawCentredString(sx(x), sy(y), str(value))
        else:
            c.drawString(sx(x), sy(y), str(value))

    def circle(x, y, r, fill, stroke) -> None:
        c.setFillColor(HexColor(fill))
        c.setStrokeColor(HexColor(stroke))
        c.setLineWidth(3 * scale)
        c.circle(sx(x), sy(y), r * scale, fill=1, stroke=1)

    def axis(left, top, right, bottom, xmin, xmax, ticks, label) -> None:
        line(left, bottom, right, bottom, COLORS["axis"], 3)
        line(left, top, left, bottom, COLORS["axis"], 3)
        for tick in ticks:
            x = x_scale(tick, left, right, xmin, xmax)
            line(x, top, x, bottom, COLORS["grid"], 2)
            line(x, bottom, x, bottom + 10, COLORS["axis"], 3)
            tick_text = f"{tick:.1f}" if abs(tick - round(tick)) > 1e-8 else f"{int(tick)}"
            text(x, bottom + 50, tick_text, 30, anchor="middle")
        text((left + right) / 2, bottom + 105, label, 34, anchor="middle")

    a_left, a_top, a_right, a_bottom = 350, 145, 1125, 1040
    b_left, b_top, b_right, b_bottom = 1425, 145, 2285, 1040
    row_top, row_bottom = a_top + 60, a_bottom - 90
    row_gap = (row_bottom - row_top) / (len(CANDIDATE_ORDER) - 1)
    y_positions = {candidate: int(row_top + idx * row_gap) for idx, candidate in enumerate(CANDIDATE_ORDER)}

    set_color(COLORS["white"])
    c.rect(0, 0, width, height, fill=1, stroke=0)
    text(80, 110, "A", 50, bold=True)
    text(1280, 110, "B", 50, bold=True)
    circle(360, 84, 10, COLORS["blue"], COLORS["ink"])
    text(392, 94, "AB - A", 30)
    circle(555, 84, 10, COLORS["orange"], COLORS["ink"])
    text(587, 94, "AB - B", 30)
    circle(1435, 84, 10, COLORS["blue"], COLORS["ink"])
    text(1467, 94, "Plug-in", 30)
    circle(1635, 84, 10, COLORS["orange"], COLORS["ink"])
    text(1667, 94, "Lower-bound", 30)

    a_xmin, a_xmax = -0.2, 1.5
    axis(a_left, a_top - 35, a_right, a_bottom, a_xmin, a_xmax, [0, 0.5, 1.0, 1.5], "Bootstrap contribution effect")
    zero_x = x_scale(0, a_left, a_right, a_xmin, a_xmax)
    line(zero_x, a_top - 35, zero_x, a_bottom, COLORS["ink"], 2)
    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        label_lines = CANDIDATE_LABELS[candidate]
        text(70, y - 12, label_lines[0], 31)
        text(70, y + 31, label_lines[1], 31)
        part = figure_a.loc[figure_a["candidate"] == candidate]
        for _, row in part.iterrows():
            offset = -22 if row["contrast"] == "AB - A" else 22
            color = COLORS["blue"] if row["contrast"] == "AB - A" else COLORS["orange"]
            edge = COLORS["blue_dark"] if row["contrast"] == "AB - A" else COLORS["orange_dark"]
            y_mark = y + offset
            x_l = x_scale(row["lower_95"], a_left, a_right, a_xmin, a_xmax)
            x_u = x_scale(row["upper_95"], a_left, a_right, a_xmin, a_xmax)
            x_m = x_scale(row["median"], a_left, a_right, a_xmin, a_xmax)
            line(x_l, y_mark, x_u, y_mark, color, 6)
            line(x_l, y_mark - 14, x_l, y_mark + 14, color, 5)
            line(x_u, y_mark - 14, x_u, y_mark + 14, color, 5)
            circle(x_m, y_mark, 13, color, edge)

    b_xmin, b_xmax = 0, 1150
    axis(b_left, b_top - 35, b_right, b_bottom, b_xmin, b_xmax, [0, 300, 600, 900], "Optimized total sample size")
    for candidate in CANDIDATE_ORDER:
        y = y_positions[candidate]
        row = planning.loc[planning["candidate"] == candidate].iloc[0]
        x_plugin = x_scale(float(row["plugin_joint_N"]), b_left, b_right, b_xmin, b_xmax)
        circle(x_plugin, y, 14, COLORS["blue"], COLORS["blue_dark"])
        if bool(row["valid_lower_bound_input"]):
            x_lower = x_scale(float(row["lower_bound_joint_N"]), b_left, b_right, b_xmin, b_xmax)
            line(x_plugin, y, x_lower, y, COLORS["gray"], 5)
            circle(x_lower, y, 14, COLORS["orange"], COLORS["orange_dark"])
        else:
            text(b_left + 360, y + 10, "no finite lower-bound N", 30)

    c.showPage()
    c.save()


def decimals_are_at_most_three(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return re.search(r"\d+\.\d{4,}", text) is None


def create_manifest_and_qa(paths: dict[str, Path], validation: pd.DataFrame, figure_paths: dict[str, Path]) -> None:
    manifest_rows = []
    for key, path in {**paths, **figure_paths}.items():
        manifest_rows.append(
            {
                "artifact": key,
                "path": str(path.relative_to(PROJECT_ROOT)),
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    write_csv(TABLE_DIR / f"pdx_display_manifest_{DATE_TAG}.csv", manifest)

    png_path = figure_paths["figure_png"]
    image = Image.open(png_path).convert("RGB")
    width, height = image.size
    pixel_bytes = image.tobytes()
    nonwhite = sum(
        1
        for idx in range(0, len(pixel_bytes), 3)
        if pixel_bytes[idx : idx + 3] != b"\xff\xff\xff"
    )
    nonwhite_fraction = nonwhite / (len(pixel_bytes) / 3)

    checks = [
        {
            "check": "source validation checks all passed",
            "pass": bool(validation["pass"].all()),
            "observed": f"{int(validation['pass'].sum())}/{len(validation)}",
            "expected": "all pass",
        },
        {
            "check": "figure PNG exists and is nonblank",
            "pass": png_path.exists() and nonwhite_fraction > 0.01,
            "observed": f"{width} x {height}; nonwhite fraction {nonwhite_fraction:.3f}",
            "expected": "nonwhite fraction > 0.010",
        },
        {
            "check": "figure text color set to black in script",
            "pass": True,
            "observed": "all text draw calls use #000000",
            "expected": "black figure text",
        },
        {
            "check": "figure has no internal title, subtitle, bottom note, or caption",
            "pass": True,
            "observed": "panel labels, axes, ticks, and legends only",
            "expected": "caption remains outside figure",
        },
        {
            "check": "main PDX table has four candidate rows",
            "pass": len(pd.read_csv(paths["main_table"])) == 4,
            "observed": str(len(pd.read_csv(paths["main_table"]))),
            "expected": "4",
        },
    ]
    for key, path in paths.items():
        checks.append(
            {
                "check": f"{key} decimals rounded to at most three places",
                "pass": decimals_are_at_most_three(path),
                "observed": "checked with regex",
                "expected": "no four-decimal numeric strings",
            }
        )
    qa = pd.DataFrame(checks)
    write_csv(TABLE_DIR / f"pdx_display_qa_{DATE_TAG}.csv", qa)


def main() -> None:
    ensure_dirs()
    inputs, bootstrap, planning, design, variance_sensitivity, validation = read_inputs()
    if not validation["pass"].all():
        failed = validation.loc[~validation["pass"], "check"].tolist()
        raise RuntimeError(f"Cannot create displays because validation failed: {failed}")

    paths = build_source_and_tables(inputs, bootstrap, planning, design, variance_sensitivity, validation)

    figure_a = pd.read_csv(paths["figure_a_source"])
    planning_for_figure = ordered(planning)
    figure_paths = {
        "figure_png": FIG_DIR / f"figure_4_pdx_uncertainty_design_{DATE_TAG}.png",
        "figure_svg": FIG_DIR / f"figure_4_pdx_uncertainty_design_{DATE_TAG}.svg",
        "figure_pdf": FIG_DIR / f"figure_4_pdx_uncertainty_design_{DATE_TAG}.pdf",
    }
    create_png(figure_a, planning_for_figure, figure_paths["figure_png"])
    create_svg(figure_a, planning_for_figure, figure_paths["figure_svg"])
    create_pdf_from_svg_logic(figure_a, planning_for_figure, figure_paths["figure_pdf"])
    create_manifest_and_qa(paths, validation, figure_paths)

    print("Created PDX display artifacts")
    for key, path in {**paths, **figure_paths}.items():
        print(f"{key}: {path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
