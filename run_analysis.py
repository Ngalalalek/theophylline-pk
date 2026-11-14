"""Theoph: NCA, individual 1-cmt fits (+/- tlag), multiple dosing.

python run_analysis.py  ->  results/*.csv, figures/*.png
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from nca import nca
from pk_models import multiple_dose, one_compartment_oral, one_compartment_oral_lag

RESULTS, FIGURES = Path("results"), Path("figures")
TARGET = (10, 20)  # mg/L
plt.rcParams["font.size"] = 9


def load_data():
    df = pd.read_csv("data/theoph.csv")
    df = df.rename(columns={"Subject": "id", "Wt": "weight", "Dose": "dose_mg_per_kg", "Time": "time"})
    df["dose_mg"] = df["dose_mg_per_kg"] * df["weight"]  # 原始数据里 Dose 是 mg/kg
    return df[["id", "weight", "dose_mg", "time", "conc"]]


def fit_subject(model, t, c, dose, p0_list, bounds):
    # sigma=c 相当于权重 1/C^2，几个初值都试一下，取 WRSS 最小的
    best = None
    for p0 in p0_list:
        try:
            p, cov = curve_fit(lambda t, *p: model(t, dose, *p), t, c, p0=p0,
                               sigma=c, bounds=bounds, maxfev=20000)
        except RuntimeError:
            continue
        wrss = np.sum(((c - model(t, dose, *p)) / c) ** 2)
        if best is None or wrss < best[2]:
            best = (p, cov, wrss)
    p, cov, wrss = best
    n, k = len(c), len(p)
    aic = n * np.log(wrss / n) + 2 * k
    se_pct = 100 * np.sqrt(np.diag(cov)) / p
    return p, se_pct, wrss, aic


def fit_all(df):
    rows = []
    for sid, g in df.groupby("id"):
        g = g[g["time"] > 0]  # 去掉给药前的点
        t, c, dose = g["time"].values, g["conc"].values, g["dose_mg"].iloc[0]

        p1, se1, wrss1, aic1 = fit_subject(one_compartment_oral, t, c, dose,
                                           [[1.5, 0.1, 30]], ([0.05, 0.01, 1], [20, 1, 200]))
        # tlag 不能超过第一个采样时间，初值都取在它下面
        p2, se2, wrss2, aic2 = fit_subject(one_compartment_oral_lag, t, c, dose,
                                           [[1.5, 0.1, 30, 0.1 * t[0]], [1.5, 0.1, 30, 0.5 * t[0]], [1.5, 0.1, 30, 0.9 * t[0]]],
                                           ([0.05, 0.01, 1, 0], [20, 1, 200, t[0]]))
        rows.append({
            "id": sid, "dose_mg": dose,
            "ka": p1[0], "ke": p1[1], "V_F": p1[2],
            "se_ka_pct": se1[0], "se_ke_pct": se1[1], "se_V_pct": se1[2],
            "wrss": wrss1, "aic": aic1,
            "lag_ka": p2[0], "lag_ke": p2[1], "lag_V_F": p2[2], "lag_tlag": p2[3],
            "lag_wrss": wrss2, "lag_aic": aic2,
        })
    fits = pd.DataFrame(rows)
    fits["CL_F"] = fits["ke"] * fits["V_F"]
    fits["t_half"] = np.log(2) / fits["ke"]
    fits["delta_aic"] = fits["lag_aic"] - fits["aic"]  # < 0: lag better
    return fits


def residuals(df, fits):
    out = []
    for sid, g in df.groupby("id"):
        g = g[g["time"] > 0]
        f = fits.set_index("id").loc[sid]
        t, c, dose = g["time"].values, g["conc"].values, f["dose_mg"]
        pred1 = one_compartment_oral(t, dose, f["ka"], f["ke"], f["V_F"])
        pred2 = one_compartment_oral_lag(t, dose, f["lag_ka"], f["lag_ke"], f["lag_V_F"], f["lag_tlag"])
        out.append(pd.DataFrame({"id": sid, "time": t, "obs": c,
                                 "res_1cmt": (c - pred1) / c, "res_lag": (c - pred2) / c}))
    return pd.concat(out, ignore_index=True)


def gmean(x):
    return float(np.exp(np.log(x).mean()))


def simulate_regimens(fits):
    # 典型患者用几何均值，ka 偏得很厉害 (subject 9 ~10/h)
    # ka, cl, V = gmean(fits["ka"]), gmean(fits["CL_F"]), gmean(fits["V_F"])
    # ka, V, tlag = gmean(fits["lag_ka"]), gmean(fits["lag_V_F"]), gmean(fits["lag_tlag"])
    ka, V = gmean(fits["lag_ka"]), gmean(fits["lag_V_F"])
    tlag = fits["lag_tlag"].median()  # 有几个人 tlag = 0，几何均值会变成 0
    cl = gmean(fits["lag_ke"] * fits["lag_V_F"])
    ke = cl / V
    # TODO 这里要不要用 tlag 模型的结果? 稳态下影响应该不大
    t = np.linspace(0, 96, 2000)
    regimens = {"300 mg every 8 h": (300, 8), "450 mg every 12 h": (450, 12)}
    curves, rows = {}, []
    for label, (dose, tau) in regimens.items():
        c = multiple_dose(t - tlag, dose, tau, int(96 / tau), ka, ke, V)
        curves[label] = c
        last = t >= 96 - tau  # 最后一个给药间隔，近似稳态
        rows.append({"regimen": label, "Css_max": c[last].max(), "Css_min": c[last].min(),
                     "Css_avg": dose / (cl * tau), "accumulation_ratio": 1 / (1 - np.exp(-ke * tau)),
                     "time_to_90pct_ss_h": np.log(10) / ke})
    return t, curves, pd.DataFrame(rows), (ka, ke, V, tlag)


def plot_observed(df):
    # 原始数据先看一下，线性和对数坐标
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    # 默认颜色只有 10 个，12 个人会重复，换 tab20
    for i, (sid, g) in enumerate(df.groupby("id")):
        col = plt.cm.tab20(i)
        axes[0].plot(g["time"], g["conc"], "o-", ms=3, lw=1, color=col, label=str(sid))
        g = g[g["conc"] > 0]  # 0 时刻浓度是 0 的点在对数坐标下画不出来
        axes[1].plot(g["time"], g["conc"], "o-", ms=3, lw=1, color=col)
    axes[1].set_yscale("log")
    for ax in axes:
        ax.set_xlabel("Time (h)")
    axes[0].set_ylabel("Conc (mg/L)")
    axes[0].legend(title="Subject", fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES / "observed_data.png", dpi=120)
    plt.close(fig)


def plot_fits(df, fits):
    fig, axes = plt.subplots(3, 4, figsize=(11, 7.5), sharex=True, sharey=True)
    tt = np.linspace(0, 25, 500)
    for ax, (sid, g) in zip(axes.flat, df.groupby("id")):
        f = fits.set_index("id").loc[sid]
        ax.plot(tt, one_compartment_oral(tt, f["dose_mg"], f["ka"], f["ke"], f["V_F"]), color="C0", label="1-cmt")
        ax.plot(tt, one_compartment_oral_lag(tt, f["dose_mg"], f["lag_ka"], f["lag_ke"], f["lag_V_F"], f["lag_tlag"]),
                color="C1", ls="--", label="1-cmt + tlag")
        ax.plot(g["time"], g["conc"], "o", color="k", ms=4, label="obs")
        ax.set_title(f"Subject {sid}")
    for ax in axes[-1]:
        ax.set_xlabel("Time (h)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Conc (mg/L)")
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURES / "individual_fits.png", dpi=120)
    plt.close(fig)


def plot_residuals(res):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, col, title in [(axes[0], "res_1cmt", "1-cmt"), (axes[1], "res_lag", "1-cmt + tlag")]:
        ax.axhline(0, color="k", lw=1, ls="--")
        ax.plot(res["time"], res[col], "o", color="C0", ms=4, alpha=0.8)
        ax.set(title=title, xlabel="Time (h)")
    axes[0].set_ylabel("(obs - pred) / obs")
    fig.tight_layout()
    fig.savefig(FIGURES / "residuals.png", dpi=120)
    plt.close(fig)


def plot_multiple_dose(t, curves):
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.axhspan(*TARGET, color="0.85", alpha=0.6, lw=0, label="10-20 mg/L")
    for (label, c), color in zip(curves.items(), ["C0", "C1"]):
        ax.plot(t, c, color=color, lw=1.6, label=label)
    ax.set(xlabel="Time (h)", ylabel="Conc (mg/L)", xlim=(0, 96), title="Typical patient, 900 mg/day")
    ax.legend(loc="lower right")
    fig.savefig(FIGURES / "multiple_dosing.png", dpi=120)
    plt.close(fig)


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    df = load_data()
    plot_observed(df)

    # NCA
    rows = []
    for sid, g in df.groupby("id"):
        r = {"id": sid}
        r.update(nca(g["time"].values, g["conc"].values, g["dose_mg"].iloc[0]))
        rows.append(r)
    nca_tab = pd.DataFrame(rows)
    nca_tab.to_csv(RESULTS / "nca_results.csv", index=False, float_format="%.3f")

    # 每个人单独拟合，1-cmt 和加 tlag 的
    fits = fit_all(df)
    print(fits.head())
    fits.to_csv(RESULTS / "model_fits.csv", index=False, float_format="%.4f")
    res = residuals(df, fits)

    # 多次给药，典型患者
    t, curves, reg, (ka, ke, V, tlag) = simulate_regimens(fits)
    reg.to_csv(RESULTS / "multiple_dosing.csv", index=False, float_format="%.2f")

    plot_fits(df, fits)
    plot_residuals(res)
    plot_multiple_dose(t, curves)

    def ms(s):
        return f"{s.mean():.3g} +/- {s.std():.2g}"
    print("NCA")
    for k in ["Cmax", "Tmax", "AUC_0_inf", "t_half", "CL_F", "V_F"]:
        print(f"  {k:10s} {ms(nca_tab[k])}")
    print("1-cmt fits")
    for k in ["ka", "ke", "V_F", "CL_F", "t_half"]:
        print(f"  {k:10s} {ms(fits[k])}   median {fits[k].median():.3g}")
    print(f"  median SE%  ka {fits['se_ka_pct'].median():.0f}  ke {fits['se_ke_pct'].median():.0f}  V {fits['se_V_pct'].median():.0f}")
    print(f"  CL/F model/NCA {(fits['CL_F'] / nca_tab['CL_F']).mean():.3f}")
    print(f"\ntypical patient (geo means): ka {ka:.2f}  ke {ke:.4f}  V {V:.1f}  CL {ke * V:.2f}  tlag {tlag:.3f}")
    print(reg.round(2).to_string(index=False))
    print()
    print(fits[["id", "lag_tlag", "aic", "lag_aic", "delta_aic"]].round(2).to_string(index=False))
    print(f"  dAIC < -2 in {(fits['delta_aic'] < -2).sum()}/12, tlag {ms(fits['lag_tlag'])} h")
    early = res["time"] < 2
    res_sum = pd.DataFrame({"model": ["1-cmt", "1-cmt + tlag"],
                            "mean_abs_res_before_2h": [res.loc[early, "res_1cmt"].abs().mean(), res.loc[early, "res_lag"].abs().mean()],
                            "mean_abs_res_after_2h": [res.loc[~early, "res_1cmt"].abs().mean(), res.loc[~early, "res_lag"].abs().mean()]})
    res_sum.to_csv(RESULTS / "residual_summary.csv", index=False, float_format="%.3f")
    print(res_sum.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
