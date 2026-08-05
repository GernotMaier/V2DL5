"""Tests for configuration merging and validation."""

from argparse import Namespace

import pytest

from v2dl5.configuration import _merge_config, configuration


def test_merge_config_preserves_nested_defaults():
    merged = _merge_config(
        {"datasets": {"safe_mask": {"methods": ["old"], "parameters": {"aeff_percent": 10}}}},
        {"datasets": {"safe_mask": {"parameters": {"aeff_percent": 5}}}},
    )

    assert merged["datasets"]["safe_mask"] == {
        "methods": ["old"],
        "parameters": {"aeff_percent": 5},
    }


def test_configuration_merges_partial_yaml(tmp_path):
    config_file = tmp_path / "config.yml"
    config_file.write_text(
        "datasets:\n  safe_mask:\n    parameters:\n      aeff_percent: 5\n",
        encoding="utf-8",
    )
    args = Namespace(config=str(config_file), run_list=None, output_dir=str(tmp_path))

    config = configuration(args)

    assert config["datasets"]["safe_mask"]["methods"] == ["aeff-default", "aeff-max"]
    assert config["datasets"]["safe_mask"]["parameters"]["aeff_percent"] == 5


def test_configuration_rejects_empty_yaml(tmp_path):
    config_file = tmp_path / "empty.yml"
    config_file.write_text("", encoding="utf-8")
    args = Namespace(config=str(config_file), run_list=None, output_dir=str(tmp_path))

    with pytest.raises(ValueError, match="empty"):
        configuration(args)
