#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
0_run_all.py
一键运行完整文献计量学流程，并将所有步骤的输出（stdout/stderr）
保存到本地日志文件，方便后续查看与复现。

日志位置：output/logs/run_YYYYMMDD_HHMMSS.log
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml


STEPS = [
    "1_parse_and_merge.py",
    "2_descriptive_analysis.py",
    "3_keyword_analysis.py",
    "4_generate_report.py",
]


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(config: dict) -> Path:
    """创建带时间戳的日志文件，返回路径。"""
    log_dir = Path(config.get("output", {}).get("root", "output")) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"run_{timestamp}.log"
    return log_path


def log_line(file_handle, line: str, to_console: bool = True):
    """同时写入日志文件并（可选）输出到控制台。"""
    file_handle.write(line)
    file_handle.flush()
    if to_console:
        print(line, end="", flush=True)


def run_command(cmd: list, log_handle, step_name: str) -> int:
    """流式运行命令，把 stdout 和 stderr 实时写入日志。"""
    log_line(log_handle, f"\n===== {step_name} =====\n")
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        encoding="utf-8",
        errors="replace",
    )
    for line in process.stdout:
        log_line(log_handle, line)
    process.wait()
    return process.returncode


def run_pipeline(config_path: str):
    config = load_config(config_path)
    log_path = setup_log_file(config)

    with open(log_path, "w", encoding="utf-8") as log_handle:
        start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = (
            f"Pipeline log started at {start_time}\n"
            f"Config: {Path(config_path).resolve()}\n"
            f"Python: {sys.executable}\n"
            f"Steps: {', '.join(STEPS)}\n"
        )
        log_line(log_handle, header)

        python = sys.executable
        for step in STEPS:
            step_path = Path(step)
            if not step_path.exists():
                msg = f"ERROR: {step} not found.\n"
                log_line(log_handle, msg)
                sys.exit(1)

            print(f"\n[0_run_all] Running {step} ...", flush=True)
            log_line(log_handle, f"\n[0_run_all] Running {step} ...\n", to_console=False)

            returncode = run_command([python, str(step_path), "--config", config_path],
                                     log_handle, step)
            if returncode != 0:
                msg = f"\n[0_run_all] {step} failed with exit code {returncode}.\n"
                log_line(log_handle, msg)
                print(f"\n详细日志见：{log_path}")
                sys.exit(returncode)

        end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        footer = (
            f"\nPipeline completed successfully at {end_time}\n"
            f"Outputs: {Path(config['output']['root']).resolve()}\n"
        )
        log_line(log_handle, footer)

    print(f"\n[0_run_all] Pipeline completed successfully.")
    print(f"[0_run_all] Log saved to: {log_path}")
    print(f"[0_run_all] Outputs are in: {Path(config['output']['root']).resolve()}")


def main():
    parser = argparse.ArgumentParser(description="Run the full bibliometric pipeline with logging.")
    parser.add_argument("--config", default="config.yaml", help="Path to config YAML")
    args = parser.parse_args()
    run_pipeline(args.config)


if __name__ == "__main__":
    main()
