package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

public class ApplyVariantAnnotationScoresIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/apply-test/test.snp.final.vcf",
                "--tranches-file", "/home/slee/working/vqsr/scalable/train-test/test.snp.tranches.csv",
                "--recal-file", "/home/slee/working/vqsr/scalable/score-test/test.snp.recal.vcf",
                "--truth-sensitivity-filter-level", "98",
                "-mode", "SNP",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/apply-test/test.indel.final.vcf",
                "--tranches-file", "/home/slee/working/vqsr/scalable/train-test/test.indel.tranches.csv",
                "--recal-file", "/home/slee/working/vqsr/scalable/score-test/test.indel.recal.vcf",
                "--truth-sensitivity-filter-level", "98",
                "-mode", "INDEL",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}