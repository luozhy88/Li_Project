#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1_parse_and_merge.py
读取 PubMed MEDLINE/CSV 与可选的 WoS CiteSpace 导出文件，
进行解析、去重，输出合并数据集与 CiteSpace 可读的 WoS 格式题录。
所有路径与参数均来自 config.yaml，脚本本身不硬编码。
"""

import argparse
import os
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ensure_dirs(config: dict):
    for d in [config["output"]["root"], config["output"]["figures"],
              config["output"]["tables"], config["output"]["report"]]:
        Path(d).mkdir(parents=True, exist_ok=True)


def parse_medline(medline_path: str) -> pd.DataFrame:
    """解析 PubMed MEDLINE 格式，返回 DataFrame。"""
    if not os.path.exists(medline_path):
        return pd.DataFrame()

    with open(medline_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    records = []
    current = None

    def flush_record(rec):
        if rec and rec.get("pmid"):
            records.append(rec)

    for raw in lines:
        line = raw.rstrip("\n")
        if line.strip() == "":
            continue
        # MEDLINE 标签为前 4 个字符，后续为内容；续行以 6 个空格开头
        if line.startswith("      ") and current is not None:
            tag = "continuation"
        else:
            tag = line[:4].strip()

        value = line[6:] if len(line) > 6 else line[5:]

        if tag == "PMID":
            flush_record(current)
            current = {"pmid": value.strip()}
        elif current is None:
            continue
        elif tag == "continuation":
            # 追加到上一个字段
            last_tag = current.get("_last_tag")
            if last_tag:
                if isinstance(current[last_tag], list):
                    current[last_tag][-1] = current[last_tag][-1] + " " + value.strip()
                else:
                    current[last_tag] = current[last_tag] + " " + value.strip()
        elif tag in ("TI", "AB", "DP", "TA", "JT", "SO", "LA", "VI",
                     "IP", "PG", "PMC", "EDAT", "MHDA", "CRDT"):
            current[tag] = value.strip()
            current["_last_tag"] = tag
        elif tag in ("AU", "FAU", "AD", "MH", "OT", "PT", "AID", "IS"):
            current.setdefault(tag, []).append(value.strip())
            current["_last_tag"] = tag
        elif tag == "VI":
            current["VI"] = value.strip()
            current["_last_tag"] = "VI"
        else:
            # 其他字段视需要保留
            current.setdefault(tag, []).append(value.strip())
            current["_last_tag"] = tag

    flush_record(current)

    if not records:
        return pd.DataFrame()

    rows = []
    for rec in records:
        pmid = rec.get("pmid", "")
        title = rec.get("TI", "")
        abstract = rec.get("AB", "")
        year = extract_year(rec.get("DP", ""))
        journal_abbr = rec.get("TA", "")
        journal_full = rec.get("JT", "")
        journal = journal_full if journal_full else journal_abbr
        authors = rec.get("AU", [])
        full_authors = rec.get("FAU", [])
        first_author = authors[0] if authors else ""
        affiliations = rec.get("AD", [])
        mesh = rec.get("MH", [])
        other_terms = rec.get("OT", [])
        publication_types = rec.get("PT", [])
        doi = extract_doi(rec.get("AID", []))
        volume = rec.get("VI", "")
        issue = rec.get("IP", "")
        issn = "; ".join(rec.get("IS", []))
        pages = rec.get("PG", "")

        rows.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "year": year,
            "journal": journal,
            "journal_abbr": journal_abbr,
            "authors": "; ".join(authors),
            "full_authors": "; ".join(full_authors),
            "first_author": first_author,
            "affiliations": " | ".join(affiliations),
            "mesh": "; ".join(mesh),
            "other_terms": "; ".join(other_terms),
            "keywords": "; ".join(clean_keywords(mesh + other_terms)),
            "publication_types": "; ".join(publication_types),
            "doi": doi,
            "volume": volume,
            "issue": issue,
            "issn": issn,
            "pages": pages,
            "source_database": "PubMed",
            "wos_id": "",
        })
    return pd.DataFrame(rows)


def extract_year(dp: str) -> int:
    if not dp:
        return None
    m = re.search(r"\b(19|20)\d{2}\b", dp)
    return int(m.group(0)) if m else None


def extract_doi(aid_list):
    for aid in aid_list:
        if "[doi]" in aid:
            return aid.replace("[doi]", "").strip()
        if aid.startswith("10."):
            return aid.strip()
    return ""


def clean_keywords(term_list):
    """清理 MeSH/OT，去掉副主题词并去重。"""
    cleaned = []
    for term in term_list:
        # MeSH 主/副主题词以 '/' 分隔，取第一部分
        root = term.split("/")[0]
        # 去掉 MeSH 主要主题词标记 *
        root = root.replace("*", "").strip()
        if root and root not in cleaned:
            cleaned.append(root)
    return cleaned


def parse_wos(wos_path: str) -> pd.DataFrame:
    """解析 WoS plain-text / CiteSpace 格式（如提供）。当前实现覆盖常见标签。"""
    if not os.path.exists(wos_path):
        return pd.DataFrame()

    with open(wos_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    records = []
    current = None

    def flush():
        if current and current.get("title"):
            records.append(current)

    for raw in lines:
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        # WoS 标签为 2 个大写字符（字母或数字，如 PT, AU, Z9, ER），后续为空格或行尾；
        # 续行通常以空格开头或无标签。
        tag_match = re.match(r"^([A-Z0-9]{2})(?:\s|$)", line)
        if tag_match:
            tag = tag_match.group(1)
            value = line[3:].strip() if len(line) > 2 else ""
        else:
            tag = "continuation"
            value = line.strip()

        if tag == "PT":
            flush()
            current = {"publication_types": value}
        elif current is None:
            continue
        elif tag == "continuation":
            last = current.get("_last_tag")
            if last:
                prev = current.get(last)
                if isinstance(prev, list):
                    # AU/AF 字段的续行是新作者，而不是上一个作者的延续
                    if last in ("authors", "full_authors"):
                        prev.append(value)
                    elif prev:
                        prev[-1] = (prev[-1] + " " + value).strip()
                elif last in ("keywords",):
                    # DE/ID 字段的续行是新关键词，用分号分隔
                    current[last] = (str(prev or "") + "; " + value).strip("; ")
                else:
                    current[last] = (str(prev or "") + " " + value).strip()
        elif tag == "AU":
            current.setdefault("authors", []).append(value)
            current["_last_tag"] = "authors"
        elif tag == "AF":
            current.setdefault("full_authors", []).append(value)
            current["_last_tag"] = "full_authors"
        elif tag == "TI":
            current["title"] = value
            current["_last_tag"] = "title"
        elif tag == "SO":
            current["journal"] = value
            current["_last_tag"] = "journal"
        elif tag == "AB":
            current["abstract"] = value
            current["_last_tag"] = "abstract"
        elif tag == "DE":
            current["keywords"] = value
            current["_last_tag"] = "keywords"
        elif tag == "DI":
            current["doi"] = value
            current["_last_tag"] = "doi"
        elif tag == "PY":
            try:
                current["year"] = int(value)
            except ValueError:
                current["year"] = None
            current["_last_tag"] = "year"
        elif tag == "VL":
            current["volume"] = value
            current["_last_tag"] = "volume"
        elif tag == "IS":
            current["issue"] = value
            current["_last_tag"] = "issue"
        elif tag == "BP":
            current["pages"] = value
            current["_last_tag"] = "pages"
        elif tag == "UT":
            current["wos_id"] = value
            current["_last_tag"] = "wos_id"
        elif tag == "C1":
            current.setdefault("affiliations", []).append(value)
            current["_last_tag"] = "affiliations"
        elif tag == "ER":
            flush()
            current = None
        else:
            # 其他 WoS 标签（LA, DT, CR, NR, TC, Z9, U1, U2, PU, PI, PA, SN, J9, PD,
            # EP, PG, WC, WE, SC, GA, OA, DA 等）：保留原始值并更新 last_tag，
            # 使后续续行不会污染上一个有效字段。
            current.setdefault(tag.lower(), []).append(value.strip())
            current["_last_tag"] = tag.lower()

    flush()

    if not records:
        return pd.DataFrame()

    rows = []
    for rec in records:
        authors = rec.get("authors", [])
        rows.append({
            "pmid": "",
            "title": rec.get("title", ""),
            "abstract": rec.get("abstract", ""),
            "year": rec.get("year"),
            "journal": rec.get("journal", ""),
            "journal_abbr": "",
            "authors": "; ".join(authors),
            "full_authors": "; ".join(rec.get("full_authors", [])),
            "first_author": authors[0] if authors else "",
            "affiliations": " | ".join(rec.get("affiliations", [])),
            "mesh": "",
            "other_terms": "",
            "keywords": rec.get("keywords", ""),
            "publication_types": rec.get("publication_types", ""),
            "doi": rec.get("doi", ""),
            "volume": rec.get("volume", ""),
            "issue": rec.get("issue", ""),
            "pages": rec.get("pages", ""),
            "source_database": "WoS",
            "wos_id": rec.get("wos_id", ""),
        })
    return pd.DataFrame(rows)


def canonicalize_author_names(df: pd.DataFrame) -> pd.DataFrame:
    """合并作者姓名拼写变体（连字符/空格/逗号/大小写差异，如 'Carrasco Moro R'
    与 'Carrasco-Moro R'、WoS 风格 'Larner, AJ' 与 PubMed 风格 'Larner AJ'）。"""
    def key(name: str) -> str:
        k = re.sub(r"[-‐‑‒–—,.'’]", " ", name.lower())
        return re.sub(r"\s+", " ", k).strip()

    variants: dict = {}
    for col in ("authors", "first_author"):
        for val in df[col].fillna(""):
            for nm in str(val).split(";"):
                nm = nm.strip()
                if nm:
                    variants.setdefault(key(nm), Counter())[nm] += 1

    display = {}
    for k, counter in variants.items():
        # 优先：出现次数多 > 含连字符 > 不带逗号 > 较短
        best = sorted(counter.items(),
                      key=lambda kv: (kv[1], "-" in kv[0], "," not in kv[0], -len(kv[0])),
                      reverse=True)[0][0]
        display[k] = best.replace(",", "")

    df = df.copy()
    for col in ("authors", "first_author"):
        df[col] = df[col].fillna("").apply(
            lambda v: "; ".join(display.get(key(nm.strip()), nm.strip())
                                for nm in str(v).split(";") if nm.strip()))
    return df


def normalize_title(title: str, stopwords: list) -> str:
    if not title:
        return ""
    text = title.lower()
    text = re.sub(r"[^a-z0-9\s]", " ", text)
    tokens = [t for t in text.split() if t and t not in stopwords]
    return " ".join(sorted(tokens))


def deduplicate(pubmed_df: pd.DataFrame, wos_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """按 PMID、DOI、归一化标题去重；PubMed 记录优先保留。"""
    stopwords = set([w.lower() for w in config["deduplication"]["title_stopwords"]])

    required_cols = ["pmid", "doi", "title", "source_database"]

    def prep(df):
        df = df.copy()
        for col in required_cols:
            if col not in df.columns:
                df[col] = ""
        df["doi"] = df["doi"].astype(str).str.lower().str.strip()
        df["pmid"] = df["pmid"].astype(str).str.strip()
        df["normalized_title"] = df["title"].apply(lambda x: normalize_title(x, stopwords))
        # 来源优先级：PubMed 0，WoS 1
        df["source_priority"] = df["source_database"].apply(lambda s: 0 if s == "PubMed" else 1)
        return df

    pubmed_df = prep(pubmed_df)
    wos_df = prep(wos_df)

    merged = pd.concat([pubmed_df, wos_df], ignore_index=True, sort=False)
    merged.sort_values(by="source_priority", inplace=True)
    merged.reset_index(drop=True, inplace=True)

    seen = {"pmid": set(), "doi": set(), "normalized_title": set()}
    keep = []
    removed = []

    for idx, row in merged.iterrows():
        dup_flag = False
        for key in config["deduplication"]["keys"]:
            val = str(row.get(key, "")).strip()
            if not val or val == "nan":
                continue
            if val in seen[key]:
                dup_flag = True
                break
        if dup_flag:
            removed.append(idx)
        else:
            keep.append(idx)
            for key in config["deduplication"]["keys"]:
                val = str(row.get(key, "")).strip()
                if val and val != "nan":
                    seen[key].add(val)

    deduped = merged.loc[keep].copy()
    deduped.drop(columns=["source_priority", "normalized_title"], inplace=True, errors="ignore")
    deduped = canonicalize_author_names(deduped)

    # 生成去重报告
    report_lines = [
        "Deduplication Report",
        "====================",
        f"PubMed records loaded: {len(pubmed_df)}",
        f"WoS records loaded:    {len(wos_df)}",
        f"Total before dedup:    {len(merged)}",
        f"Records removed:       {len(removed)}",
        f"Unique records kept:   {len(deduped)}",
        "",
        "Source distribution after deduplication:",
        deduped["source_database"].value_counts().to_string(),
    ]
    report_path = Path(config["output"]["root"]) / "01_deduplication_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print(f"[1_parse_and_merge] Deduplication report saved to {report_path}")
    return deduped


def export_citespace(df: pd.DataFrame, output_path: str):
    """将合并结果导出为 CiteSpace 可读的 WoS plain-text 格式。"""
    lines = []
    for _, row in df.iterrows():
        pt = "J" if "Journal Article" in str(row.get("publication_types", "")) else "J"
        au_list = [a.strip() for a in str(row.get("authors", "")).split(";") if a.strip()]
        af_list = [a.strip() for a in str(row.get("full_authors", "")).split(";") if a.strip()]
        title = str(row.get("title", "")).replace("\n", " ")
        journal = str(row.get("journal", "")).replace("\n", " ")
        abstract = str(row.get("abstract", "")).replace("\n", " ")
        keywords = str(row.get("keywords", "")).replace(";", "; ")
        affiliations = str(row.get("affiliations", "")).replace("\n", " ").replace(" | ", "; ")
        year = ""
        if pd.notna(row.get("year")):
            try:
                year = str(int(float(row["year"])))
            except (ValueError, TypeError):
                year = str(row["year"])
        doi = str(row.get("doi", ""))
        volume = str(row.get("volume", ""))
        issue = str(row.get("issue", ""))
        pages = str(row.get("pages", ""))
        wos_id = str(row.get("wos_id", "")) or (f"PMID:{row.get('pmid')}" if row.get("pmid") else "")

        lines.append(f"PT {pt}")
        for au in au_list:
            lines.append(f"AU {au}")
        for af in af_list:
            lines.append(f"AF {af}")
        lines.append(f"TI {title}")
        lines.append(f"SO {journal}")
        lines.append(f"LA English")
        lines.append(f"DT Article")
        if keywords.strip():
            lines.append(f"DE {keywords}")
        if abstract.strip():
            lines.append(f"AB {abstract}")
        if affiliations.strip():
            lines.append(f"C1 {affiliations}")
        lines.append(f"PY {year}")
        if volume:
            lines.append(f"VL {volume}")
        if issue:
            lines.append(f"IS {issue}")
        if pages:
            lines.append(f"BP {pages}")
        if doi:
            lines.append(f"DI {doi}")
        lines.append(f"UT {wos_id}")
        lines.append("ER")
        lines.append("")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[1_parse_and_merge] CiteSpace input saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Parse PubMed/WoS and merge/deduplicate.")
    parser.add_argument("--config", default="config.yaml", help="Path to config YAML")
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_dirs(config)

    # 读取 PubMed
    print("[1_parse_and_merge] Parsing PubMed MEDLINE...")
    pubmed_df = parse_medline(config["input"]["pubmed_medline"])
    if pubmed_df.empty:
        print("[1_parse_and_merge] WARNING: No PubMed records parsed.")
    else:
        print(f"[1_parse_and_merge] Parsed {len(pubmed_df)} PubMed records.")

    # 读取可选 WoS
    print("[1_parse_and_merge] Checking WoS file...")
    wos_df = parse_wos(config["input"]["wos_file"])
    if wos_df.empty:
        print(f"[1_parse_and_merge] WoS file not found or empty: {config['input']['wos_file']}")
        print("[1_parse_and_merge] Continuing with PubMed records only.")
    else:
        print(f"[1_parse_and_merge] Parsed {len(wos_df)} WoS records.")

    # 去重
    merged = deduplicate(pubmed_df, wos_df, config)

    # 保存合并数据
    merged_csv = Path(config["output"]["root"]) / "01_merged_records.csv"
    merged.to_csv(merged_csv, index=False, encoding="utf-8-sig")
    print(f"[1_parse_and_merge] Merged records saved to {merged_csv}")

    # 导出 CiteSpace 格式
    export_citespace(merged, config["citespace"]["output_file"])


if __name__ == "__main__":
    main()
