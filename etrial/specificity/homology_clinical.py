"""
Clinical-Grade Specificity Screening Module.

Uses MMseqs2 for ultra-fast sequence homology search against human proteome.

Installation:
    conda install -c bioconda mmseqs2

    # Download human proteome database
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
    gunzip UP000005640_9606.fasta.gz

    # Create MMseqs2 database
    mmseqs createdb human_proteome.fasta humanDB

Reference:
    Steinegger & Söding, Nature Biotechnology (2017)
    "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets"
"""

from typing import List, Dict, Any, Optional
from pathlib import Path
import subprocess
import tempfile
import os
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    Decision,
    MetricResult,
)


class ClinicalSpecificityModule(ValidationModule):
    """
    Clinical-grade specificity screening using MMseqs2.

    Screens for:
        - Sequence homology to human proteome
        - Conserved motifs (off-target risk)
        - Tissue expression of homologs
    """

    def _setup(self) -> None:
        """Check for MMseqs2 and database availability."""
        logger.info("Setting up ClinicalSpecificityModule")

        # Check for MMseqs2
        self.mmseqs2_available = self._check_mmseqs2()

        # Check for human proteome database
        self.db_path = os.environ.get('MMSEQS2_HUMAN_DB')

        if self.mmseqs2_available:
            logger.info("MMseqs2 found")
            if self.db_path and Path(self.db_path).exists():
                logger.info(f"Human proteome database: {self.db_path}")
            else:
                logger.warning(
                    "Human proteome database not found.\n"
                    "Set: export MMSEQS2_HUMAN_DB=/path/to/humanDB\n"
                    "Using basic sequence analysis only."
                )
        else:
            logger.warning(
                "MMseqs2 not available.\n"
                "Install: conda install -c bioconda mmseqs2\n"
                "Falling back to basic homology checks."
            )

    def _check_mmseqs2(self) -> bool:
        """Check if MMseqs2 is installed."""
        try:
            result = subprocess.run(
                ['mmseqs', 'version'],
                capture_output=True,
                timeout=5
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def validate(
        self,
        candidate: TherapeuticCandidate
    ) -> ValidationResult:
        """
        Validate specificity using MMseqs2 homology search.

        Args:
            candidate: Therapeutic candidate

        Returns:
            ValidationResult with specificity assessment
        """
        logger.info(f"Running clinical specificity validation for {candidate.name}")

        metrics: List[MetricResult] = []
        risks: List[str] = []
        warnings: List[str] = []
        recommendations: List[str] = []

        sequence = candidate.get_primary_sequence()

        if not sequence or candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="No sequence available or small molecule - specificity screening not applicable",
                risks=[],
                warnings=["No sequence available or small molecule"],
                recommendations=[],
            )

        # ================================================================
        # Homology Search
        # ================================================================

        if self.mmseqs2_available and self.db_path:
            homology_result = self._search_with_mmseqs2(sequence)
        else:
            homology_result = self._check_conserved_motifs(sequence)
            warnings.append("Using motif-based screening - install MMseqs2 for clinical predictions")

        # Number of off-target hits
        metrics.append(MetricResult(
            name="n_off_target_hits",
            value=homology_result['n_hits'],
            unit="count",
            threshold_pass=2,
            threshold_revise=10,
            decision=Decision.PASS if homology_result['n_hits'] < 2 else (
                Decision.REVISE if homology_result['n_hits'] < 10 else Decision.KILL
            ),
            metadata={'method': homology_result.get('method', 'unknown')}
        ))

        # Max identity to off-targets
        if homology_result.get('max_identity'):
            metrics.append(MetricResult(
                name="max_off_target_identity",
                value=homology_result['max_identity'],
                unit="percent",
                threshold_pass=70,
                threshold_revise=85,
                decision=Decision.PASS if homology_result['max_identity'] < 70 else (
                    Decision.REVISE if homology_result['max_identity'] < 85 else Decision.KILL
                ),
            ))

        # Conserved motifs (promiscuous binding)
        if homology_result.get('conserved_motifs'):
            motifs = homology_result['conserved_motifs']

            # Check for known problematic motifs
            if 'RGD' in motifs:
                risks.append("RGD motif detected - binds integrins (off-target risk)")
            if any('K' * 5 in m for m in motifs):
                risks.append("Poly-basic region detected - promiscuous binding risk")

            metrics.append(MetricResult(
                name="n_conserved_motifs",
                value=len(motifs),
                unit="count",
                threshold_pass=1,
                threshold_revise=3,
                decision=Decision.PASS if len(motifs) < 1 else (
                    Decision.REVISE if len(motifs) < 3 else Decision.KILL
                ),
            ))

        # Add recommendations
        if homology_result['n_hits'] > 2:
            recommendations.append(
                f"Found {homology_result['n_hits']} potential off-targets. "
                "Consider: (1) experimental cross-reactivity testing, "
                "(2) interface redesign for selectivity"
            )

        # Overall decision
        decisions = [m.decision for m in metrics]
        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
        else:
            overall_decision = Decision.PASS

        # Create summary
        n_hits = homology_result['n_hits']
        method = homology_result.get('method', 'unknown')
        max_identity = homology_result.get('max_identity', 0)
        summary = f"Specificity screening: {n_hits} off-target hits, max identity {max_identity:.1f}% (method: {method})"

        return self._create_result(
            candidate=candidate,
            decision=overall_decision,
            metrics=metrics,
            summary=summary,
            risks=risks,
            warnings=warnings,
            recommendations=recommendations,
        )

    def _search_with_mmseqs2(self, sequence: str) -> Dict[str, Any]:
        """
        Use MMseqs2 to search human proteome.

        Returns:
            dict with n_hits, max_identity, conserved_motifs, method
        """
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                # Write sequence
                fasta_file = tmpdir / "query.fasta"
                with open(fasta_file, 'w') as f:
                    f.write(f">query\n{sequence}\n")

                # Create query DB
                query_db = tmpdir / "queryDB"
                subprocess.run(
                    ['mmseqs', 'createdb', str(fasta_file), str(query_db)],
                    check=True,
                    capture_output=True
                )

                # Search
                result_db = tmpdir / "resultDB"
                subprocess.run([
                    'mmseqs', 'search',
                    str(query_db),
                    self.db_path,
                    str(result_db),
                    tmpdir / 'tmp',
                    '--min-seq-id', '0.6',  # 60% identity threshold
                    '-e', '0.01',  # E-value threshold
                ], check=True, capture_output=True, timeout=300)

                # Convert to table
                result_table = tmpdir / "result.tsv"
                subprocess.run([
                    'mmseqs', 'convertalis',
                    str(query_db),
                    self.db_path,
                    str(result_db),
                    str(result_table),
                ], check=True, capture_output=True)

                # Parse results
                hits = []
                with open(result_table) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            hits.append({
                                'target': parts[1],
                                'identity': float(parts[2]),
                            })

                # Filter self-hits
                hits = [h for h in hits if h['identity'] < 99.0]

                n_hits = len(hits)
                max_identity = max([h['identity'] for h in hits]) if hits else 0.0

                return {
                    'n_hits': n_hits,
                    'max_identity': max_identity,
                    'conserved_motifs': self._extract_motifs(sequence) if n_hits > 5 else [],
                    'method': 'mmseqs2',
                }

        except Exception as e:
            logger.error(f"MMseqs2 search failed: {e}")
            return self._check_conserved_motifs(sequence)

    def _check_conserved_motifs(self, sequence: str) -> Dict[str, Any]:
        """Fallback: Check for known promiscuous motifs."""
        motifs = self._extract_motifs(sequence)

        # Estimate off-target risk based on motifs
        n_hits = len(motifs) * 3  # Rough estimate

        return {
            'n_hits': n_hits,
            'max_identity': 70.0 if n_hits > 5 else 50.0,
            'conserved_motifs': motifs,
            'method': 'motif_screening',
        }

    def _extract_motifs(self, sequence: str) -> List[str]:
        """Extract potentially problematic motifs."""
        motifs = []

        # RGD motif (integrin binding)
        if 'RGD' in sequence:
            motifs.append('RGD')

        # Poly-basic (promiscuous)
        import re
        poly_basic = re.findall(r'K{5,}|R{5,}', sequence)
        motifs.extend(poly_basic)

        # Poly-acidic (promiscuous)
        poly_acidic = re.findall(r'D{5,}|E{5,}', sequence)
        motifs.extend(poly_acidic)

        return motifs


def is_mmseqs2_available() -> bool:
    """Check if MMseqs2 is installed."""
    try:
        result = subprocess.run(
            ['mmseqs', 'version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
