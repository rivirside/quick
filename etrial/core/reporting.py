"""
Report generation for eTrial validation results.

Generates publication-quality dossiers in multiple formats (PDF, HTML, JSON)
documenting validation results, decisions, and recommendations.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from jinja2 import Template, Environment, FileSystemLoader
from loguru import logger

from etrial.core.base import PipelineResult, ValidationResult, Decision


class ReportGenerator:
    """
    Generates validation dossiers in multiple formats.

    Produces comprehensive reports including:
    - Executive summary
    - Module-by-module results
    - Visualizations
    - Decision recommendations
    - Risk mitigation strategies
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize report generator.

        Args:
            config: Reporting configuration dictionary
        """
        self.config = config or {}
        self.template_dir = Path(__file__).parent.parent / "templates"
        self.template_dir.mkdir(parents=True, exist_ok=True)

        # Set up Jinja2 environment
        self.jinja_env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=True,
        )

    def generate_report(
        self,
        pipeline_result: PipelineResult,
        output_path: Path,
        format: str = "html",
    ) -> Path:
        """
        Generate validation dossier.

        Args:
            pipeline_result: Pipeline result to document
            output_path: Path for output file
            format: Output format ('pdf', 'html', 'json')

        Returns:
            Path to generated report
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Generating {format.upper()} report: {output_path}")

        if format == "json":
            return self._generate_json(pipeline_result, output_path)
        elif format == "html":
            return self._generate_html(pipeline_result, output_path)
        elif format == "pdf":
            return self._generate_pdf(pipeline_result, output_path)
        else:
            raise ValueError(f"Unsupported format: {format}")

    def _generate_json(
        self,
        pipeline_result: PipelineResult,
        output_path: Path
    ) -> Path:
        """Generate JSON report."""
        json_data = pipeline_result.to_dict()

        with open(output_path, 'w') as f:
            json.dump(json_data, f, indent=2, default=str)

        logger.info(f"JSON report saved: {output_path}")
        return output_path

    def _generate_html(
        self,
        pipeline_result: PipelineResult,
        output_path: Path
    ) -> Path:
        """Generate HTML report."""
        # Create HTML template if it doesn't exist
        template_path = self.template_dir / "report_template.html"
        if not template_path.exists():
            self._create_default_html_template(template_path)

        template = self.jinja_env.get_template("report_template.html")

        # Prepare data for template
        context = self._prepare_report_context(pipeline_result)

        # Render template
        html_content = template.render(**context)

        # Save HTML
        with open(output_path, 'w') as f:
            f.write(html_content)

        logger.info(f"HTML report saved: {output_path}")
        return output_path

    def _generate_pdf(
        self,
        pipeline_result: PipelineResult,
        output_path: Path
    ) -> Path:
        """Generate PDF report (via HTML -> PDF)."""
        # First generate HTML
        html_path = output_path.with_suffix('.html')
        self._generate_html(pipeline_result, html_path)

        try:
            # Try using weasyprint for PDF generation
            from weasyprint import HTML
            HTML(filename=str(html_path)).write_pdf(output_path)
            html_path.unlink()  # Remove intermediate HTML
            logger.info(f"PDF report saved: {output_path}")

        except ImportError:
            logger.warning(
                "weasyprint not installed, PDF generation unavailable. "
                "HTML report saved instead."
            )
            # Rename HTML to PDF extension (it's still HTML)
            html_path.rename(output_path.with_suffix('.html'))

        return output_path

    def _prepare_report_context(self, pipeline_result: PipelineResult) -> Dict[str, Any]:
        """Prepare data context for report template."""
        candidate = pipeline_result.candidate

        # Group modules by decision
        modules_pass = []
        modules_revise = []
        modules_kill = []

        for module_name, result in pipeline_result.module_results.items():
            module_summary = {
                "name": module_name,
                "decision": result.decision.value,
                "summary": result.summary,
                "metrics": result.metrics,
                "recommendations": result.recommendations,
                "risks": result.risks,
                "runtime": result.runtime_seconds,
            }

            if result.decision == Decision.PASS:
                modules_pass.append(module_summary)
            elif result.decision == Decision.REVISE:
                modules_revise.append(module_summary)
            elif result.decision == Decision.KILL:
                modules_kill.append(module_summary)

        # Collect all recommendations and risks
        all_recommendations = []
        all_risks = []

        for result in pipeline_result.module_results.values():
            all_recommendations.extend(result.recommendations)
            all_risks.extend(result.risks)

        # Executive summary statistics
        exec_summary = {
            "overall_decision": pipeline_result.overall_decision.value,
            "candidate_name": candidate.name,
            "modality": candidate.modality.value,
            "target": candidate.target,
            "total_modules": len(pipeline_result.module_results),
            "modules_pass": len(modules_pass),
            "modules_revise": len(modules_revise),
            "modules_kill": len(modules_kill),
            "total_runtime": pipeline_result.total_runtime_seconds,
            "timestamp": pipeline_result.timestamp.strftime("%Y-%m-%d %H:%M:%S"),
        }

        context = {
            "exec_summary": exec_summary,
            "candidate": candidate,
            "modules_pass": modules_pass,
            "modules_revise": modules_revise,
            "modules_kill": modules_kill,
            "all_recommendations": all_recommendations,
            "all_risks": all_risks,
            "pipeline_result": pipeline_result,
        }

        return context

    def _create_default_html_template(self, template_path: Path) -> None:
        """Create default HTML template."""
        html_template = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>eTrial Validation Dossier - {{ exec_summary.candidate_name }}</title>
    <style>
        body {
            font-family: 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 8px;
            margin-bottom: 30px;
        }
        h1 { margin: 0; font-size: 2.5em; }
        h2 { color: #667eea; border-bottom: 2px solid #667eea; padding-bottom: 10px; }
        h3 { color: #764ba2; }

        .decision-badge {
            display: inline-block;
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 1.2em;
            margin: 10px 0;
        }
        .decision-PASS { background-color: #10b981; color: white; }
        .decision-REVISE { background-color: #f59e0b; color: white; }
        .decision-KILL { background-color: #ef4444; color: white; }
        .decision-INFORMATIVE { background-color: #6b7280; color: white; }

        .section {
            background: white;
            padding: 25px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .exec-summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .stat-card {
            background: #f9fafb;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #667eea;
        }
        .stat-label { font-size: 0.9em; color: #6b7280; }
        .stat-value { font-size: 1.8em; font-weight: bold; color: #111827; }

        .module-result {
            border-left: 4px solid #e5e7eb;
            padding-left: 20px;
            margin: 15px 0;
        }
        .module-result.PASS { border-left-color: #10b981; }
        .module-result.REVISE { border-left-color: #f59e0b; }
        .module-result.KILL { border-left-color: #ef4444; }

        .metric-table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        .metric-table th {
            background-color: #f3f4f6;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            border-bottom: 2px solid #d1d5db;
        }
        .metric-table td {
            padding: 10px 12px;
            border-bottom: 1px solid #e5e7eb;
        }

        .recommendation-list, .risk-list {
            list-style-type: none;
            padding: 0;
        }
        .recommendation-list li, .risk-list li {
            padding: 10px;
            margin: 8px 0;
            border-radius: 6px;
        }
        .recommendation-list li {
            background-color: #dbeafe;
            border-left: 4px solid #3b82f6;
        }
        .risk-list li {
            background-color: #fee2e2;
            border-left: 4px solid #ef4444;
        }

        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #6b7280;
            font-size: 0.9em;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>eTrial Validation Dossier</h1>
        <p style="font-size: 1.2em; margin: 10px 0 0 0;">
            {{ exec_summary.candidate_name }} - {{ exec_summary.modality }}
        </p>
    </div>

    <!-- Executive Summary -->
    <div class="section">
        <h2>Executive Summary</h2>

        <p><strong>Overall Decision:</strong>
            <span class="decision-badge decision-{{ exec_summary.overall_decision }}">
                {{ exec_summary.overall_decision }}
            </span>
        </p>

        <div class="exec-summary">
            <div class="stat-card">
                <div class="stat-label">Target</div>
                <div class="stat-value">{{ exec_summary.target }}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Modules Run</div>
                <div class="stat-value">{{ exec_summary.total_modules }}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Passed</div>
                <div class="stat-value" style="color: #10b981;">{{ exec_summary.modules_pass }}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Revise</div>
                <div class="stat-value" style="color: #f59e0b;">{{ exec_summary.modules_revise }}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Killed</div>
                <div class="stat-value" style="color: #ef4444;">{{ exec_summary.modules_kill }}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Runtime</div>
                <div class="stat-value">{{ "%.1f"|format(exec_summary.total_runtime) }}s</div>
            </div>
        </div>

        <p><strong>Timestamp:</strong> {{ exec_summary.timestamp }}</p>
    </div>

    <!-- Modules: KILL -->
    {% if modules_kill %}
    <div class="section">
        <h2 style="color: #ef4444;">Critical Issues (KILL)</h2>
        {% for module in modules_kill %}
        <div class="module-result KILL">
            <h3>{{ module.name }}</h3>
            <p><strong>Decision:</strong> <span class="decision-badge decision-KILL">KILL</span></p>
            <p>{{ module.summary }}</p>

            {% if module.risks %}
            <h4>Risks:</h4>
            <ul class="risk-list">
                {% for risk in module.risks %}
                <li>{{ risk }}</li>
                {% endfor %}
            </ul>
            {% endif %}
        </div>
        {% endfor %}
    </div>
    {% endif %}

    <!-- Modules: REVISE -->
    {% if modules_revise %}
    <div class="section">
        <h2 style="color: #f59e0b;">Addressable Issues (REVISE)</h2>
        {% for module in modules_revise %}
        <div class="module-result REVISE">
            <h3>{{ module.name }}</h3>
            <p><strong>Decision:</strong> <span class="decision-badge decision-REVISE">REVISE</span></p>
            <p>{{ module.summary }}</p>

            {% if module.recommendations %}
            <h4>Recommendations:</h4>
            <ul class="recommendation-list">
                {% for rec in module.recommendations %}
                <li>{{ rec }}</li>
                {% endfor %}
            </ul>
            {% endif %}
        </div>
        {% endfor %}
    </div>
    {% endif %}

    <!-- Modules: PASS -->
    {% if modules_pass %}
    <div class="section">
        <h2 style="color: #10b981;">Passed Modules</h2>
        {% for module in modules_pass %}
        <div class="module-result PASS">
            <h3>{{ module.name }}</h3>
            <p><strong>Decision:</strong> <span class="decision-badge decision-PASS">PASS</span></p>
            <p>{{ module.summary }}</p>
        </div>
        {% endfor %}
    </div>
    {% endif %}

    <!-- All Recommendations -->
    {% if all_recommendations %}
    <div class="section">
        <h2>All Recommendations</h2>
        <ul class="recommendation-list">
            {% for rec in all_recommendations %}
            <li>{{ rec }}</li>
            {% endfor %}
        </ul>
    </div>
    {% endif %}

    <div class="footer">
        <p>Generated by eTrial v0.1.0 | {{ exec_summary.timestamp }}</p>
        <p>In-Silico Clinical Trial Platform</p>
    </div>
</body>
</html>
"""

        with open(template_path, 'w') as f:
            f.write(html_template)

        logger.info(f"Default HTML template created: {template_path}")

    def generate_executive_summary(
        self,
        pipeline_result: PipelineResult
    ) -> str:
        """
        Generate 1-page executive summary text.

        Args:
            pipeline_result: Pipeline validation result

        Returns:
            Executive summary text
        """
        candidate = pipeline_result.candidate
        decision = pipeline_result.overall_decision

        summary_parts = [
            f"EXECUTIVE SUMMARY",
            f"=" * 60,
            f"",
            f"Candidate: {candidate.name}",
            f"Modality: {candidate.modality.value}",
            f"Target: {candidate.target}",
            f"",
            f"OVERALL DECISION: {decision.value}",
            f"",
        ]

        # Module summary
        pass_count = len([r for r in pipeline_result.module_results.values()
                         if r.decision == Decision.PASS])
        revise_count = len([r for r in pipeline_result.module_results.values()
                           if r.decision == Decision.REVISE])
        kill_count = len([r for r in pipeline_result.module_results.values()
                         if r.decision == Decision.KILL])

        summary_parts.extend([
            f"Modules Run: {len(pipeline_result.module_results)}",
            f"  - PASS: {pass_count}",
            f"  - REVISE: {revise_count}",
            f"  - KILL: {kill_count}",
            f"",
            f"Runtime: {pipeline_result.total_runtime_seconds:.2f} seconds",
            f"Timestamp: {pipeline_result.timestamp.strftime('%Y-%m-%d %H:%M:%S')}",
            f"",
        ])

        # Key findings
        if kill_count > 0:
            summary_parts.append("CRITICAL ISSUES (KILL):")
            for name, result in pipeline_result.module_results.items():
                if result.decision == Decision.KILL:
                    summary_parts.append(f"  - {name}: {result.summary}")
            summary_parts.append("")

        if revise_count > 0:
            summary_parts.append("ADDRESSABLE ISSUES (REVISE):")
            for name, result in pipeline_result.module_results.items():
                if result.decision == Decision.REVISE:
                    summary_parts.append(f"  - {name}: {result.summary}")
            summary_parts.append("")

        # Recommendations
        all_recs = []
        for result in pipeline_result.module_results.values():
            all_recs.extend(result.recommendations)

        if all_recs:
            summary_parts.append("TOP RECOMMENDATIONS:")
            for i, rec in enumerate(all_recs[:5], 1):
                summary_parts.append(f"  {i}. {rec}")

        return "\n".join(summary_parts)
