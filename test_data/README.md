# Test Data Directory

This directory is intended for storing test MGF (Mascot Generic Format) files to validate the pipeline functionality.

## Adding Test Data

1. Place your `.mgf` test files in this directory
2. Run the pipeline with:
   ```bash
   nextflow run ../main.nf --input 'test_data/*.mgf' -profile test,docker
   ```

## MGF File Format

MGF files should contain mass spectrometry data in the following format:

```
BEGIN IONS
TITLE=Sample_Spectrum_1
RTINSECONDS=123.45
PEPMASS=500.25 100.0
CHARGE=2+
123.4567 1000.0
234.5678 2000.0
345.6789 1500.0
END IONS

BEGIN IONS
TITLE=Sample_Spectrum_2
...
END IONS
```

## Expected Output

For each `.mgf` file processed, you should see:
- Annotation TSV file in `results/annotations/`
- HTML QC report in `results/qc_reports/`

## Sample Test Data

If you need sample MGF files for testing:

1. **Option A**: Use publicly available datasets
   - GNPS (Global Natural Products Social Molecular Networking)
   - MassIVE (Mass Spectrometry Interactive Virtual Environment)
   - MetaboLights

2. **Option B**: Generate minimal test file manually

Create a minimal test file `test_sample.mgf`:

```mgf
BEGIN IONS
TITLE=Test_Spectrum_1
RTINSECONDS=60.0
PEPMASS=180.0634 1000.0
CHARGE=1+
59.0133 100.0
89.0239 200.0
180.0634 1000.0
END IONS

BEGIN IONS
TITLE=Test_Spectrum_2
RTINSECONDS=120.0
PEPMASS=256.2402 800.0
CHARGE=1+
112.1120 150.0
140.1070 300.0
256.2402 800.0
END IONS
```

## Troubleshooting

If the pipeline fails with test data:
1. Verify MGF file format is correct
2. Check that files have `.mgf` extension
3. Ensure files are not empty
4. Validate spectrum data is properly formatted

## Notes

- Larger test files will take longer to process
- Start with small test files (2-10 spectra) for quick validation
- Test data is excluded from Docker builds (see `.dockerignore`)
