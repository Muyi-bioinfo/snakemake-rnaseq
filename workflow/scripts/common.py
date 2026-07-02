#!/usr/bin/env python3
"""
RNA-seq Snakemake 流程 - Python 辅助函数
支持双端(paired)和单端(single)两种读段模式，通过 config.yaml 切换
所有输出路径统一从 config["paths"]["outputs"] 读取，改名只需改 config
"""

import csv
import json
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Any, Union


# ============================================================
# 内部辅助
# ============================================================

def _out(config: dict, step: str) -> str:
    """快捷取值: _out(config, "fastqc_raw") → "results/01_fastqc_raw" """
    return config["paths"]["outputs"][step]


def _get_read_mode(config: dict) -> str:
    return config.get("read_pattern", {}).get("mode", "paired")


# ============================================================
# 样本表加载（CSV/TSV/TXT，根据扩展名自动选择分隔符）
# ============================================================

def load_samples_table(filepath: str) -> Tuple[List[str], Dict[str, List[str]]]:
    """从 CSV/TSV/TXT 文件加载样本表，根据扩展名自动选择分隔符

    文件必须包含 header 行，至少含 ``sample`` 和 ``group`` 两列。
    返回 ``(samples_list, groups_dict)``。

    >>> load_samples_table("config/samples.csv")
    (["WT_1", "WT_2", "KO_1"], {"WT": ["WT_1", "WT_2"], "KO": ["KO_1"]})
    """
    ext = Path(filepath).suffix.lower()
    delimiter = "," if ext == ".csv" else "\t"

    samples: List[str] = []
    groups: Dict[str, List[str]] = {}
    with open(filepath, newline="") as f:
        # DictReader 自动将第一行解析为列名
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            s = row["sample"].strip()
            g = row["group"].strip()
            samples.append(s)
            groups.setdefault(g, []).append(s)
    return samples, groups


# ============================================================
# 基础工具函数
# ============================================================

def get_samples(config: dict) -> List[str]:
    """获取样本列表 — 优先从 samples_file 加载，否则回退到 YAML 内 samples"""
    samples_file = config.get("samples_file", "")
    if samples_file:
        samples, _ = load_samples_table(samples_file)
        return samples
    return config.get("samples", [])


def get_read_pairs(
    sample: str, raw_dir: str,
    r1_suffix: str = "_R1.fastq.gz", r2_suffix: str = "_R2.fastq.gz",
) -> Tuple[str, str]:
    r1 = os.path.join(raw_dir, f"{sample}{r1_suffix}")
    r2 = os.path.join(raw_dir, f"{sample}{r2_suffix}")
    return r1, r2


def to_json_str(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


# ============================================================
# 模式感知的输入/输出函数
# ============================================================

def get_read_input_list(wildcards, config: dict) -> List[str]:
    mode = _get_read_mode(config)
    raw_dir = config["paths"]["raw_dir"]
    pattern = config["read_pattern"]
    if mode == "single":
        return [os.path.join(raw_dir, f"{wildcards.sample}{pattern['single_suffix']}")]
    r1, r2 = get_read_pairs(wildcards.sample, raw_dir, pattern["r1_suffix"], pattern["r2_suffix"])
    return [r1, r2]


def get_read_inputs(wildcards, config: dict) -> dict:
    mode = _get_read_mode(config)
    raw_dir = config["paths"]["raw_dir"]
    pattern = config["read_pattern"]
    if mode == "single":
        return {"read": os.path.join(raw_dir, f"{wildcards.sample}{pattern['single_suffix']}")}
    r1, r2 = get_read_pairs(wildcards.sample, raw_dir, pattern["r1_suffix"], pattern["r2_suffix"])
    return {"r1": r1, "r2": r2}


def get_trimmed_reads(wildcards, config: dict) -> dict:
    mode = _get_read_mode(config)
    tdir = _out(config, "trimmomatic")
    if mode == "single":
        return {"read": os.path.join(tdir, f"{wildcards.sample}.trimmed.fastq.gz")}
    return {
        "r1": os.path.join(tdir, f"{wildcards.sample}_R1.trimmed.fastq.gz"),
        "r2": os.path.join(tdir, f"{wildcards.sample}_R2.trimmed.fastq.gz"),
    }


def get_trimmed_input_list(wildcards, config: dict) -> List[str]:
    mode = _get_read_mode(config)
    tdir = _out(config, "trimmomatic")
    if mode == "single":
        return [os.path.join(tdir, f"{wildcards.sample}.trimmed.fastq.gz")]
    return [
        os.path.join(tdir, f"{wildcards.sample}_R1.trimmed.fastq.gz"),
        os.path.join(tdir, f"{wildcards.sample}_R2.trimmed.fastq.gz"),
    ]


def get_hisat2_input(wildcards, config: dict) -> dict:
    mode = _get_read_mode(config)
    tdir = _out(config, "trimmomatic")
    if mode == "single":
        return {"read": os.path.join(tdir, f"{wildcards.sample}.trimmed.fastq.gz")}
    return {
        "r1": os.path.join(tdir, f"{wildcards.sample}_R1.trimmed.fastq.gz"),
        "r2": os.path.join(tdir, f"{wildcards.sample}_R2.trimmed.fastq.gz"),
    }


def get_bam(wildcards, config: dict) -> str:
    return os.path.join(_out(config, "hisat2"), f"{wildcards.sample}.sorted.bam")


# ============================================================
# FastQC 输出列表
# ============================================================

def get_fastqc_raw_outputs(wildcards, config: dict) -> List[str]:
    out = _out(config, "fastqc_raw")
    if _get_read_mode(config) == "single":
        return [f"{out}/{wildcards.sample}_fastqc.html", f"{out}/{wildcards.sample}_fastqc.zip"]
    return [
        f"{out}/{wildcards.sample}_R1_fastqc.html", f"{out}/{wildcards.sample}_R2_fastqc.html",
        f"{out}/{wildcards.sample}_R1_fastqc.zip",  f"{out}/{wildcards.sample}_R2_fastqc.zip",
    ]


def get_fastqc_trimmed_outputs(wildcards, config: dict) -> List[str]:
    out = _out(config, "fastqc_trimmed")
    if _get_read_mode(config) == "single":
        return [f"{out}/{wildcards.sample}.trimmed_fastqc.html", f"{out}/{wildcards.sample}.trimmed_fastqc.zip"]
    return [
        f"{out}/{wildcards.sample}_R1.trimmed_fastqc.html", f"{out}/{wildcards.sample}_R2.trimmed_fastqc.html",
        f"{out}/{wildcards.sample}_R1.trimmed_fastqc.zip",  f"{out}/{wildcards.sample}_R2.trimmed_fastqc.zip",
    ]


# ============================================================
# Trimmomatic 输出列表
# ============================================================

def get_trimmomatic_outputs(wildcards, config: dict) -> List[str]:
    out = _out(config, "trimmomatic")
    if _get_read_mode(config) == "single":
        return [f"{out}/{wildcards.sample}.trimmed.fastq.gz"]
    return [
        f"{out}/{wildcards.sample}_R1.trimmed.fastq.gz",
        f"{out}/{wildcards.sample}_R2.trimmed.fastq.gz",
        f"{out}/{wildcards.sample}_R1.unpaired.fastq.gz",
        f"{out}/{wildcards.sample}_R2.unpaired.fastq.gz",
    ]


# ============================================================
# 顶层 target 汇总（供 rule all 使用）
# ============================================================

def get_all_pipeline_targets(config: dict, samples: List[str]) -> List[str]:
    targets = []
    mode = _get_read_mode(config)
    O = config["paths"]["outputs"]               # shorthand

    for s in samples:
        if mode == "single":
            targets.append(f"{O['fastqc_raw']}/{s}_fastqc.html")
        else:
            targets.extend([
                f"{O['fastqc_raw']}/{s}_R1_fastqc.html",
                f"{O['fastqc_raw']}/{s}_R2_fastqc.html",
            ])

    for s in samples:
        if mode == "single":
            targets.append(f"{O['fastqc_trimmed']}/{s}.trimmed_fastqc.html")
        else:
            targets.extend([
                f"{O['fastqc_trimmed']}/{s}_R1.trimmed_fastqc.html",
                f"{O['fastqc_trimmed']}/{s}_R2.trimmed_fastqc.html",
            ])

    for s in samples:
        targets.append(f"{O['hisat2']}/{s}.sorted.bam")

    targets.append(f"{O['featurecounts']}/featurecounts.txt")
    targets.append(f"{O['multiqc']}/multiqc_report.html")
    targets.extend(get_deseq2_outputs(config))
    targets.extend(get_clusterprofiler_outputs(config))

    return targets


# ============================================================
# DESeq2 辅助函数
# ============================================================

def get_deseq2_groups_json(config: dict) -> str:
    """获取分组 JSON — 优先从 samples_file 加载，否则回退到 YAML 内 groups"""
    samples_file = config.get("samples_file", "")
    if samples_file:
        _, groups = load_samples_table(samples_file)
        return to_json_str(groups)
    return to_json_str(config.get("groups", {}))


def get_deseq2_contrasts_json(config: dict) -> str:
    return to_json_str(config.get("deseq2", {}).get("contrasts", []))


def get_deseq2_outputs(config: dict) -> List[str]:
    deseq2_cfg = config.get("deseq2", {})
    contrasts = deseq2_cfg.get("contrasts", [])
    outdir = _out(config, "deseq2")
    outputs = [
        f"{outdir}/PCA_plot.pdf", f"{outdir}/sample_distance_heatmap.pdf",
        f"{outdir}/DEG_heatmap.pdf",
    ]
    for ct in contrasts:
        n = ct["name"]
        outputs.extend([
            f"{outdir}/{n}_all_results.csv",  f"{outdir}/{n}_significant.csv",
            f"{outdir}/{n}_MA_plot.pdf",      f"{outdir}/{n}_volcano_plot.pdf",
        ])
    return outputs


# ============================================================
# clusterProfiler 辅助函数
# ============================================================

def get_clusterprofiler_outputs(config: dict) -> List[str]:
    """返回 clusterProfiler 预期输出文件列表（按对比展开）"""
    outdir = _out(config, "clusterprofiler")
    contrasts = config.get("deseq2", {}).get("contrasts", [])
    outputs = []
    for ct in contrasts:
        n = ct["name"]
        outputs.extend([
            f"{outdir}/{n}_GO_enrichment.csv",
            f"{outdir}/{n}_GO_dotplot.pdf",
            f"{outdir}/{n}_KEGG_enrichment.csv",
            f"{outdir}/{n}_KEGG_dotplot.pdf",
        ])
    return outputs


# ============================================================
# 配置校验（启动前检查关键字段，避免跑到一半才报拼写错误）
# ============================================================

def validate_config(config: dict) -> None:
    """校验 config 关键字段的合法性，发现问题立即退出（退出码 1）"""
    errors: List[str] = []

    # 1. read_pattern.mode 必须是 paired 或 single
    mode = config.get("read_pattern", {}).get("mode", "paired")
    if mode not in ("paired", "single"):
        errors.append(f"read_pattern.mode 必须为 'paired' 或 'single'，当前值: '{mode}'")

    # 2. 模式对应后缀必须存在
    rp = config.get("read_pattern", {})
    if mode == "paired":
        if "r1_suffix" not in rp:
            errors.append("双端模式下 read_pattern 缺少 r1_suffix")
        if "r2_suffix" not in rp:
            errors.append("双端模式下 read_pattern 缺少 r2_suffix")
    if mode == "single":
        if "single_suffix" not in rp:
            errors.append("单端模式下 read_pattern 缺少 single_suffix")

    # 3. 参考文件存在性
    ref = config.get("reference", {})
    for key, label in [
        ("genome_fasta", "参考基因组 FASTA"),
        ("annotation_gtf", "注释 GTF"),
    ]:
        if key in ref:
            p = Path(ref[key])
            if not p.is_file():
                errors.append(f"reference.{key} 文件不存在: {ref[key]}")

    # 4. deseq2 contrasts 结构
    contrasts = config.get("deseq2", {}).get("contrasts", [])
    if not contrasts:
        errors.append("deseq2.contrasts 为空，至少需要定义一个对比")
    for i, ct in enumerate(contrasts):
        for field in ("name", "case", "control"):
            if field not in ct:
                errors.append(f"deseq2.contrasts[{i}] 缺少 '{field}' 字段")
        if "case" in ct and "control" in ct and ct["case"] == ct["control"]:
            errors.append(
                f"contrasts[{i}] '{ct['name']}': case 和 control 相同 ('{ct['control']}')"
            )

    # 5. clusterprofiler 段（如果存在）
    cp = config.get("clusterprofiler", {})
    if cp:
        if "org_db" not in cp:
            errors.append("clusterprofiler.org_db 未设置（如 org.Hs.eg.db）")
        if "kegg_organism" not in cp:
            errors.append("clusterprofiler.kegg_organism 未设置（如 hsa）")

    # 6. 样本表文件（如果使用 samples_file）
    samples_file = config.get("samples_file", "")
    if samples_file and not Path(samples_file).is_file():
        errors.append(f"samples_file 指向的文件不存在: {samples_file}")

    # 汇总报错
    if errors:
        print("ERROR: config.yaml 配置校验失败：", file=sys.stderr)
        for e in errors:
            print(f"  ✗ {e}", file=sys.stderr)
        raise SystemExit(1)

    print(f"✓ 配置校验通过（{mode} 模式，{len(contrasts)} 个对比）", file=sys.stderr)


# ============================================================
# 输入文件校验（启动前检查，避免跑到一半才发现文件缺失）
# ============================================================

def validate_input_files(config: dict, samples: List[str]) -> None:
    """验证所有样本的原始 FASTQ 文件是否存在。

    缺失时打印清晰的错误信息并退出（退出码 1），
    全部存在时打印确认信息。
    """
    raw_dir = Path(config["paths"]["raw_dir"])
    mode = config.get("read_pattern", {}).get("mode", "paired")
    pattern = config["read_pattern"]
    missing: List[str] = []

    for s in samples:
        if mode == "single":
            f = raw_dir / f"{s}{pattern['single_suffix']}"
            if not f.is_file():
                missing.append(str(f))
        else:
            r1 = raw_dir / f"{s}{pattern['r1_suffix']}"
            r2 = raw_dir / f"{s}{pattern['r2_suffix']}"
            if not r1.is_file():
                missing.append(str(r1))
            if not r2.is_file():
                missing.append(str(r2))

    if missing:
        print("ERROR: 以下原始 FASTQ 文件不存在：", file=sys.stderr)
        for f in sorted(missing):
            print(f"  - {f}", file=sys.stderr)
        print(f"\n请检查 {config.get('samples_file', 'config/samples.tsv')} "
              f"中的样本名是否与实际文件匹配", file=sys.stderr)
        raise SystemExit(1)

    print(f"✓ 已确认 {len(samples)} 个样本的输入文件全部存在", file=sys.stderr)
