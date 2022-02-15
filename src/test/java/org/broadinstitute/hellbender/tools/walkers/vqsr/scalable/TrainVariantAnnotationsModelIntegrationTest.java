package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

public class TrainVariantAnnotationsModelIntegrationTest extends CommandLineProgramTest {

    private static final String PYTHON_SCRIPT = packageMainResourcesDir + "tools/walkers/vqsr/scalable/isolation-forest.py";

    @Test
    public void test1kgp50ExomesCombined() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.indel",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPAS() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.as.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.snp.as",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMAll() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm.all",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMCombined() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm.all",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--ignore-variant-type",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMIndel() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.train",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
                "--mode", "SNP",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}