# CancerDataServices-CatchERR
This script will take a CDS metadata manifest file and try to blindly fix the most common errors before the validation step.

At this time, the script will catch errors for case in enumerated values, white space in values, and incomplete urls for files.

To run the script on a CDS template, run the following command in a terminal where R is installed for help.

```
Rscript --vanilla CDS-CatchERR.R --help
```

```
Usage: CDS-CatchERR.R [options]

CDS-CatchERR v2.0.7

Options:
	-f CHARACTER, --file=CHARACTER
		dataset file (.xlsx, .tsv, .csv)

	-t CHARACTER, --template=CHARACTER
		dataset template file, CDS_submission_metadata_template.xlsx

	-h, --help
		Show this help message and exit
```

To run this script on an example file, please use the following:

```
Rscript --vanilla CDS-CatchERR.R -t test_files_v1.3/c_fail_missing_values-v1.3.1.xlsx -f test_files_v1.3/c_fail_missing_values-v1.3.1.xlsx


Process Complete.

The output file can be found here: CancerDataServices-CatchERR/test_files_v1.3/

```

This will return an output file log that will tell you what changes were made, as well as a new manifest that will have the new changes applied.

```
ERROR: gender property contains a value that is not recognized: male
	The value in gender was changed: male ---> Male
ERROR: gender property contains a value that is not recognized: female
	The value in gender was changed: female ---> Female
PASS: race property contains all valid values.
PASS: library_source property contains all valid values.
ERROR: platform property contains a value that is not recognized: ilLUMINA
	The value in platform was changed: ilLUMINA ---> ILLUMINA
ERROR: platform property contains a value that is not recognized: ILLUMina
	The value in platform was changed: ILLUMina ---> ILLUMINA
PASS: sample_type property contains all valid values.
PASS: disease_type property contains all valid values.
```

An example for an incomplete url can be seen in the following code:

```
Rscript --vanilla CDS-CatchERR.R -t test_files_v1.3/c_fail_missing_values-v1.3.1.xlsx -f test_files_v1.3/d_fail_url_file_missing-v1.3.1.xlsx 


Process Complete.

The output file can be found here: CancerDataServices-CatchERR/test_files_v1.3/
```

With the output looking like this:

```
PASS: gender property contains all valid values.
PASS: race property contains all valid values.
PASS: study_data_types property contains all valid values.
PASS: library_source property contains all valid values.
PASS: platform property contains all valid values.
PASS: sample_type property contains all valid values.

WARNING: The file location for the file, file2.bam, has been changed: s3://bucket/ ---> s3://bucket/file2.bam
WARNING: The file location for the file, file3.bam, has been changed: s3://bucket ---> s3://bucket/file3.bam
```
