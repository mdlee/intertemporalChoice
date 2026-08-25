#!/usr/bin/env python3
"""Write overnight-mcmc-live.canvas.tsx with live Stan sequential detail."""

from __future__ import annotations

import datetime as dt
import re
import sys
from pathlib import Path

MODELS = Path(__file__).resolve().parent
DEFAULT_CANVAS = Path(
    "/home/mdlee/.cursor/projects/"
    "home-mdlee-Dropbox-GitHub-intertemporalChoice/canvases/"
    "overnight-mcmc-live.canvas.tsx"
)

FAMILIES = [
    ("exponentialExecutionHierarchical_entrop", "Exponential"),
    ("hyperbolicExecutionHierarchical_entrop", "Hyperbolic"),
    ("proportionalDifferencesExecutionHierarchical_entrop", "Prop. differences"),
    ("directDifferencesExecutionHierarchical_entrop", "Direct differences"),
    ("tradeoffExecutionHierarchical_entrop", "Tradeoff"),
    ("unifiedTradeoffExecutionHierarchical_entrop", "Unified tradeoff"),
    ("itchExecutionHierarchical_entrop", "ITCH"),
    ("hyperboloidExecutionHierarchical_entrop", "Hyperboloid"),
]


def mu_prec_stem(base: str, lab: str) -> str:
    return re.sub(r"_entrop$", f"MuPrec{lab}_entrop", base)


def stan_jobs() -> list[tuple[str, str, str]]:
    jobs = [(stem, fam, "base") for stem, fam in FAMILIES]
    for stem, fam in FAMILIES:
        jobs.append((mu_prec_stem(stem, "Half"), fam, "μ-prec ×½"))
        jobs.append((mu_prec_stem(stem, "Double"), fam, "μ-prec ×2"))
    return jobs


MIX_JOBS = [
    ("latentMixtureHierarchicalPrecision_entrop", "Latent mixture", "base"),
    ("latentMixtureHierarchicalPrecisionMuPrecHalf_entrop", "Latent mixture", "μ-prec ×½"),
    ("latentMixtureHierarchicalPrecisionMuPrecDouble_entrop", "Latent mixture", "μ-prec ×2"),
]


def classify_mix(name, st_msg, mix_msg, hier_done, has_ck, idx, mix_idx):
    met = parse_metrics(st_msg)
    if not hier_done:
        return "queued", met
    return classify_stan(name, st_msg, mix_msg, has_ck, idx, mix_idx)


def read_status(path: Path) -> tuple[str, str]:
    if not path.is_file():
        return "", ""
    lines = [ln.strip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    if not lines:
        return "", ""
    return lines[0], " ".join(lines[1:])


def first_num(msg: str, keys: list[str]):
    for key in keys:
        m = re.search(re.escape(key) + r"([0-9.Inf]+)", msg)
        if not m:
            continue
        raw = m.group(1)
        if raw.lower() == "inf":
            return float("inf")
        try:
            return float(raw)
        except ValueError:
            continue
    return None


def parse_metrics(msg: str) -> dict:
    keep_m = re.search(r"keep=([0-9]+/[0-9]+)", msg)
    return {
        "rhat": first_num(msg, ["R-hat=", "R-hat max=", "trailing R-hat=", "subset R-hat max="]),
        "best": first_num(msg, ["best R-hat=", "best="]),
        "stacked": first_num(msg, ["stacked=", "Stacked="]),
        "batch": first_num(msg, ["batch="]),
        "trail": first_num(msg, ["trail=", "trailing n="]),
        "keep": keep_m.group(1) if keep_m else "",
    }


def parse_job_index(msg: str) -> int | None:
    m = re.search(r"\[(\d+)/(\d+)\]", msg)
    return int(m.group(1)) if m else None


def parse_stan_chains(tmp_dir: Path) -> list[dict]:
    if not tmp_dir.is_dir():
        return []
    files = sorted(
        tmp_dir.glob("samples_*.csv"),
        key=lambda p: int(re.search(r"(\d+)", p.stem).group(1) or 0),
    )
    chains = []
    for path in files:
        text = path.read_text(errors="replace")
        rec = {
            "id": int(re.search(r"(\d+)", path.stem).group(1) or 0),
            "phase": "waiting",
            "draws": 0,
            "warmup": 1000,
            "samples": 1000,
            "started": "—",
            "bytes": path.stat().st_size,
        }
        w = re.search(r"num_warmup\s*=\s*(\d+)", text)
        s = re.search(r"num_samples\s*=\s*(\d+)", text)
        st = re.search(r"start_datetime\s*=\s*([^\n\r]+)", text)
        if w:
            rec["warmup"] = int(w.group(1))
        if s:
            rec["samples"] = int(s.group(1))
        if st:
            rec["started"] = st.group(1).strip()
        n_data = 0
        seen_hdr = False
        for ln in text.splitlines():
            if not ln or ln.startswith("#"):
                continue
            if not seen_hdr and "lp__" in ln:
                seen_hdr = True
                continue
            if seen_hdr:
                n_data += 1
        rec["draws"] = n_data
        if rec["samples"] > 0 and n_data >= rec["samples"]:
            rec["phase"] = "done"
        elif n_data > 0:
            rec["phase"] = "sampling"
        else:
            rec["phase"] = "warmup"
        chains.append(rec)
    return chains


def finite(x) -> bool:
    return x is not None and x == x and x != float("inf")


def fmt_r(x) -> str:
    return f"{x:.3f}" if finite(x) else "—"


def fmt_n(x) -> str:
    if not finite(x):
        return "—"
    return str(int(x)) if float(x).is_integer() else f"{x:g}"


def rhat_tone(x):
    if not finite(x):
        return None
    if x <= 1.05:
        return "success"
    if x <= 1.15:
        return "warning"
    return "danger"


def status_tone(status: str):
    return {
        "done": "success",
        "skip": "success",
        "running": "warning",
        "sampling": "warning",
        "warmup": "warning",
        "paused": "info",
        "checkpoint": "danger",
        "stopped": "danger",
        "queued": "info",
    }.get(status)


def jsx_text(x: str) -> str:
    s = str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return re.sub(r"\r?\n", " ", s)


def jsx_attr(x: str) -> str:
    return jsx_text(x).replace('"', "&quot;")


def js_str(x: str) -> str:
    return str(x).replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")


def recent_log(path: Path, n: int, skip: list[str]) -> str:
    if not path.is_file():
        return ""
    lines = [ln.strip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    keep = [ln for ln in lines if not any(s in ln for s in skip)]
    if not keep:
        return ""
    chunk = " · ".join(keep[-n:])
    return chunk[-900:]


def classify_stan(name, st_msg, seq_msg, has_ck, idx, cur_idx) -> tuple[str, dict]:
    met = parse_metrics(st_msg)
    low_st = st_msg.lower()
    low_seq = seq_msg.lower()
    seq_hit = name in seq_msg
    if "finished" in low_seq or "using saved" in low_st or "converged=1" in st_msg:
        return "done", met
    if "stop file" in low_seq and seq_hit:
        return "stopped", met
    fitting = seq_hit and ("fitting" in low_seq or "running" in low_seq)
    if fitting or (cur_idx is not None and idx == cur_idx):
        return "running", met
    if cur_idx is not None and idx < cur_idx:
        return "done", met
    if has_ck or (finite(met["stacked"]) and met["stacked"] > 0):
        return "checkpoint", met
    return "pending", met


def jsx_stat(value: str, label: str, tone=None) -> str:
    if tone:
        return f'        <Stat value="{jsx_attr(value)}" label="{jsx_attr(label)}" tone="{tone}" />'
    return f'        <Stat value="{jsx_attr(value)}" label="{jsx_attr(label)}" />'


def row_tone_js(status: str) -> str:
    t = status_tone(status)
    return "undefined" if t is None else f'"{t}"'


def write_canvas(canvas_path: Path) -> None:
    jobs = stan_jobs()
    jags_log = MODELS / "logs"
    stan_log = MODELS / "stan" / "logs"
    stan_tmp = MODELS / "stan" / "tmp"
    _, seq_j = read_status(jags_log / "runHierarchicalExecutionSequential.status")
    _, mix = read_status(jags_log / "runLatentMixtureRobustnessSequential.status")
    _, seq_s = read_status(stan_log / "runHierarchicalExecutionSequentialStan.status")
    _, mix_s = read_status(stan_log / "runLatentMixtureSequentialStan.status")
    hier_done = "finished at" in seq_s.lower()
    mix_idx = parse_job_index(mix_s) if hier_done else None

    cur_idx = parse_job_index(seq_s)
    n_stan = len(jobs)
    rows = []
    tones = []
    n_done = 0
    cur_name = ""
    live_name = ""
    cur_status = ""
    cur_met = parse_metrics("")

    for j, (name, family, variant) in enumerate(jobs, start=1):
        _, st_msg = read_status(stan_log / f"{name}.status")
        has_ck = (stan_log / f"{name}.checkpoint.mat").is_file()
        status, met = classify_stan(name, st_msg, seq_s, has_ck, j, cur_idx)
        if status in {"done", "skip"}:
            n_done += 1
        if status in {"running", "warmup", "sampling"}:
            cur_name = f"{family} ({variant})"
            live_name = name
            cur_status = status
            cur_met = met
        cells = [
            f'"{j}"',
            f'"{js_str(family)}"',
            f'"{js_str(variant)}"',
            f'"{js_str(status)}"',
            f'"{fmt_n(met["batch"])}"',
            f'"{fmt_n(met["stacked"])}"',
            f'"{fmt_n(met["trail"])}"',
            f'"{fmt_r(met["rhat"])}"',
            f'"{fmt_r(met["best"])}"',
            f'"{js_str(met["keep"] or "—")}"',
        ]
        rows.append("    [" + ", ".join(cells) + "]")
        tones.append(row_tone_js(status))

    if not live_name and cur_idx is not None and 1 <= cur_idx <= n_stan:
        live_name, fam, var = jobs[cur_idx - 1]
        cur_name = f"{fam} ({var})"
        _, st_msg = read_status(stan_log / f"{live_name}.status")
        cur_met = parse_metrics(st_msg)
        if not cur_status:
            cur_status = "running"

    mix_rows = []
    mix_tones = []
    for j, (name, family, variant) in enumerate(MIX_JOBS, start=1):
        _, st_msg = read_status(stan_log / f"{name}.status")
        has_ck = (stan_log / f"{name}.checkpoint.mat").is_file()
        status, met = classify_mix(name, st_msg, mix_s, hier_done, has_ck, j, mix_idx)
        if status in {"running", "warmup", "sampling"}:
            live_name = name
            cur_name = f"{family} ({variant})"
            cur_status = status
            cur_met = met
        mix_rows.append(
            "    ["
            + ", ".join(
                [
                    f'"{j}"',
                    f'"{js_str(family)}"',
                    f'"{js_str(variant)}"',
                    f'"{js_str(status)}"',
                    f'"{fmt_n(met["batch"])}"',
                    f'"{fmt_n(met["stacked"])}"',
                    f'"{fmt_n(met["trail"])}"',
                    f'"{fmt_r(met["rhat"])}"',
                    f'"{fmt_r(met["best"])}"',
                    f'"{js_str(met["keep"] or "—")}"',
                ]
            )
            + "]"
        )
        mix_tones.append(row_tone_js(status))

    chains = parse_stan_chains(stan_tmp / live_name) if live_name else []
    n_chains = len(chains)
    n_warm, n_samp = 1000, 1000
    draws_min = draws_max = None
    phase = "waiting"
    if chains:
        draws = [c["draws"] for c in chains]
        draws_min, draws_max = min(draws), max(draws)
        n_warm, n_samp = chains[0]["warmup"], chains[0]["samples"]
        phases = [c["phase"] for c in chains]
        if "sampling" in phases:
            phase = "sampling"
        elif "warmup" in phases:
            phase = "warmup"
        elif all(p == "done" for p in phases):
            phase = "batch complete"
        else:
            phase = phases[0]
        if cur_status in {"running", ""}:
            cur_status = phase

    jags_stopped = any(w in seq_j.lower() for w in ("stopped", "halted"))
    jags_done = "finished" in seq_j.lower()
    mix_running = bool(mix_s) and "finished at" not in mix_s.lower() and hier_done
    if "finished at" in mix_s.lower():
        stage, callout = "Stan mixture sequential — finished", "success"
    elif mix_running:
        stage, callout = f"Stan mixture sequential — {phase}", "warning"
    elif "finished" in seq_s.lower():
        stage, callout = "Stan sequential — finished (mixture queued)", "info"
    elif "stop file" in seq_s.lower():
        stage, callout = "Stan sequential — stopped", "warning"
    elif seq_s:
        stage, callout = f"Stan sequential — {phase}", "warning"
    else:
        stage, callout = "Stan sequential", "info"

    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    title = f"{stage} · Trinity CmdStan · {now}"
    body = (
        "JAGS hierarchical sequential is stopped. Stan is the live overnight queue: "
        "24 hierarchical jobs (8 families × base, μ-prec ×½, μ-prec ×2), then the "
        "11-component latent mixture (8 hierarchical cognitive models + Guess/LL/SS) "
        "and its μ-precision ×½ / ×2 robustness jobs. Engine Trinity/CmdStan, "
        f"6 chains, warmup {n_warm:g} then collect {n_samp:g} per batch, thin 1. "
        "Hierarchical gate: subset max R-hat below 1.05 and n at least 10,000 on a best 4-of-6 subset. "
        "Mixture gate: max participant z R-hat at most 1.10 on a best 4-of-6 subset, n at least 10,000. "
        "CmdStan save_warmup is false, so CSV draws stay at 0 through warmup and grow only in sampling. "
        f"Now: {cur_name or 'queue'} — {mix_s or seq_s or 'starting'}."
    )

    n_run = 1 if (cur_status in {"running", "warmup", "sampling"} or "fitting" in seq_s.lower() or "fitting" in mix_s.lower()) else 0
    if "stop file" in seq_s.lower() or "stop file" in mix_s.lower():
        n_run = 0
    if "finished at" in mix_s.lower():
        n_run = 0
    elif "finished" in seq_s.lower() and not mix_running:
        n_run = 0
    n_run = min(n_run, max(0, n_stan - n_done)) if not mix_running else n_run
    n_pend = max(0, n_stan - n_done - n_run)

    stan_recent = recent_log(
        stan_log / "runHierarchicalExecutionSequentialStan.log",
        10,
        ["No canonical to archive"],
    )
    mix_recent = recent_log(
        stan_log / "runLatentMixtureSequentialStan.log",
        8,
        ["No canonical to archive"],
    )
    run_recent = recent_log(
        jags_log / "mcmcRuns.status.log",
        8,
        ["tradeoffExecutionHierarchical_entrop"],
    )

    L: list[str] = []
    a = L.append
    a("import {")
    a("  Callout,")
    a("  Card,")
    a("  CardBody,")
    a("  CardHeader,")
    a("  H1,")
    a("  H2,")
    a("  H3,")
    a("  Row,")
    a("  Stack,")
    a("  Stat,")
    a("  Table,")
    a("  Text,")
    a("  UsageBar,")
    a('} from "cursor/canvas";')
    a("")
    a("export default function OvernightMcmcLive() {")
    a("  return (")
    a("    <Stack gap={24}>")
    a("      <Stack gap={8}>")
    a("        <H1>Overnight Stan sequential</H1>")
    a(f'        <Text tone="secondary">{jsx_text(title)}</Text>')
    a("      </Stack>")
    a(f'      <Callout tone="{callout}" title="{jsx_attr(stage)}">{jsx_text(body)}</Callout>')
    a('      <Row gap={24} align="end">')
    a(jsx_stat(cur_name or "—", "Now running"))
    a(jsx_stat(phase, "CmdStan phase"))
    a(jsx_stat(fmt_r(cur_met["rhat"]), "Subset R-hat", rhat_tone(cur_met["rhat"])))
    a(jsx_stat(fmt_r(cur_met["best"]), "Best R-hat"))
    a(jsx_stat(fmt_n(cur_met["stacked"]), "Stacked / chain"))
    a(jsx_stat(fmt_n(cur_met["batch"]), "Batch"))
    a(jsx_stat(fmt_n(cur_met["trail"]), "Trail window"))
    a(jsx_stat(cur_met["keep"] or "—", "Keep chains"))
    a(jsx_stat(f"{n_done} / {n_stan}", "Stan jobs done"))
    a("      </Row>")

    segs = []
    if n_done:
        segs.append(f'        {{ id: "done", value: {n_done}, color: "green" }}')
    if n_run:
        segs.append(f'        {{ id: "run", value: {n_run}, color: "orange" }}')
    if n_pend:
        segs.append(f'        {{ id: "pend", value: {n_pend}, color: "gray" }}')
    a(
        f'      <UsageBar total={{{n_stan}}} topLeftLabel="Stan queue" '
        f'topRightLabel="{n_done} done · {n_run} in flight · {n_pend} pending" segments={{['
    )
    a(",\n".join(segs))
    a("      ]} />")

    a("      <H2>Sequential status</H2>")
    a(f"      <Text>{jsx_text('Stan hierarchical: ' + (seq_s or '(none)'))}</Text>")
    a(f"      <Text>{jsx_text('Stan mixture: ' + (mix_s or ('queued after hierarchical' if not hier_done else '(none)')))}</Text>")
    jags_line = seq_j or ("stopped" if jags_stopped else "(idle)")
    if jags_done:
        jags_line = seq_j
    a(f'      <Text tone="secondary">{jsx_text("JAGS: " + jags_line)}</Text>')
    a(
        f'      <Text tone="secondary">{jsx_text("JAGS mixture robustness: " + (mix or "not started (JAGS queue halted)"))}</Text>'
    )

    if chains:
        chain_rows = []
        chain_tones = []
        for ch in chains:
            chain_rows.append(
                "    ["
                + ", ".join(
                    [
                        f'"{ch["id"]}"',
                        f'"{js_str(ch["phase"])}"',
                        f'"{fmt_n(ch["draws"])}"',
                        f'"{fmt_n(ch["warmup"])}"',
                        f'"{fmt_n(ch["samples"])}"',
                        f'"{js_str(ch["started"])}"',
                        f'"{fmt_n(ch["bytes"])}"',
                    ]
                )
                + "]"
            )
            chain_tones.append(
                {"done": '"success"', "sampling": '"warning"', "warmup": '"info"'}.get(
                    ch["phase"], "undefined"
                )
            )
        a("      <Card>")
        a(
            f'        <CardHeader trailing="{jsx_attr(f"{n_chains} chains · draws {fmt_n(draws_min)}-{fmt_n(draws_max)}")}">'
            f"{jsx_text(live_name or 'CmdStan batch')}</CardHeader>"
        )
        a("        <CardBody>")
        a("          <Stack gap={12}>")
        a('            <Text tone="secondary" size="small">')
        a("              Live CmdStan CSV under models/stan/tmp. Header-only files mean warmup;")
        a("              sampling appends one row per saved draw (warmup is not saved).")
        a("            </Text>")
        a("            <Table")
        a(
            '              headers={["Chain", "Phase", "Saved draws", "Warmup", "Sample", "Started (UTC)", "CSV bytes"]}'
        )
        a('              columnAlign={["right", "left", "right", "right", "right", "left", "right"]}')
        a(f'              rowTone={{[{", ".join(chain_tones)}]}}')
        a("              striped")
        a("              rows={[")
        a(",\n".join(chain_rows))
        a("              ]}")
        a("            />")
        a("          </Stack>")
        a("        </CardBody>")
        a("      </Card>")

    job_label = str(cur_idx) if cur_idx is not None else "—"
    a("      <H2>Stan job board</H2>")
    a('      <Text tone="secondary" size="small">')
    a(
        f"        {jsx_text(f'Job {job_label} of {n_stan}.')} Same 8 families as JAGS, then μ-precision ×½ and ×2. Mixture jobs start after these 24 finish."
    )
    a("      </Text>")
    a("      <Table")
    a(
        '        headers={["#", "Family", "Variant", "Status", "Batch", "Stacked", "Trail", "R-hat", "Best", "Keep"]}'
    )
    a(
        '        columnAlign={["right", "left", "left", "left", "right", "right", "right", "right", "right", "right"]}'
    )
    a(f'        rowTone={{[{", ".join(tones)}]}}')
    a("        striped")
    a("        stickyHeader")
    a("        rows={[")
    a(",\n".join(rows))
    a("        ]}")
    a("      />")

    mix_label = str(mix_idx) if mix_idx is not None else "—"
    a("      <H2>Stan latent mixture</H2>")
    a('      <Text tone="secondary" size="small">')
    a(
        f"        {jsx_text(f'Mixture job {mix_label} of {len(MIX_JOBS)}.')} "
        "Manuscript 11-component mixture (participant-level z, 8 hierarchical cognitive models + Guess/LL/SS), then μ-precision ×½ and ×2. Queued until hierarchical Stan sequential finishes."
    )
    a("      </Text>")
    a("      <Table")
    a(
        '        headers={["#", "Family", "Variant", "Status", "Batch", "Stacked", "Trail", "R-hat", "Best", "Keep"]}'
    )
    a(
        '        columnAlign={["right", "left", "left", "left", "right", "right", "right", "right", "right", "right"]}'
    )
    a(f'        rowTone={{[{", ".join(mix_tones)}]}}')
    a("        striped")
    a("        rows={[")
    a(",\n".join(mix_rows))
    a("        ]}")
    a("      />")

    if stan_recent:
        a("      <H2>Stan sequential log</H2>")
        a(f'      <Text tone="secondary" size="small">{jsx_text(stan_recent)}</Text>')
    if mix_recent:
        a("      <H2>Stan mixture sequential log</H2>")
        a(f'      <Text tone="secondary" size="small">{jsx_text(mix_recent)}</Text>')
    if run_recent:
        a("      <H3>Shared MCMC status log</H3>")
        a(f'      <Text tone="secondary" size="small">{jsx_text(run_recent)}</Text>')

    a("    </Stack>")
    a("  );")
    a("}")
    a("")

    canvas_path.parent.mkdir(parents=True, exist_ok=True)
    canvas_path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CANVAS
    write_canvas(dest)
